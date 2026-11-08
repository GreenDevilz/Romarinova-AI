"""
Dock all rosemary compounds against the prepared receptor.
Uses AutoDock Vina via subprocess and Meeko for ligand preparation.
"""
import json
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "models" / "docking_config.json"
RECEPTOR_PATH = ROOT / "data" / "3L4W_receptor.pdbqt"
CATALOGUE_PATH = ROOT / "data" / "cultivar_compound_catalogue.csv"
OUTPUT_PATH = ROOT / "results" / "docking_results.csv"

# --- Vina Executable Path ---
VINA_EXE = ROOT / "vina.exe"

# --- Temporary Directory for PDBQT files ---
TEMP_DIR = ROOT / "data" / "temp_docking"
TEMP_DIR.mkdir(parents=True, exist_ok=True)

def prepare_ligand(smiles, out_pdbqt_path):
    """Generate 3D conformer and write PDBQT to file."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return False
    mol = Chem.AddHs(mol)
    
    # Locked prep: ETKDGv3, random seed 42, MMFF94 optimize
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        return False
    AllChem.MMFFOptimizeMolecule(mol)
    
    preparer = MoleculePreparation()
    mol_setups = preparer.prepare(mol)
    pdbqt_string, is_ok, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
    
    if is_ok:
        out_pdbqt_path.write_text(pdbqt_string)
        return True
    return False

def parse_vina_affinity(out_pdbqt_path):
    """Parse the top pose affinity from the Vina output PDBQT."""
    if not out_pdbqt_path.exists():
        return None
    with open(out_pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT:"):
                # Format: REMARK VINA RESULT:    -6.2      0.000      0.000
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        return float(parts[3])
                    except ValueError:
                        pass
    return None

def main():
    if not VINA_EXE.exists():
        raise FileNotFoundError(
            f"Could not find {VINA_EXE}. "
            "Please make sure vina.exe is in the root of your project directory."
        )

    if not RECEPTOR_PATH.exists() or RECEPTOR_PATH.stat().st_size == 0:
        raise FileNotFoundError(
            f"Receptor file {RECEPTOR_PATH} is missing or empty. "
            "Did you run scripts/prepare_receptor.py?"
        )

    config = json.loads(CONFIG_PATH.read_text())
    df = pd.read_csv(CATALOGUE_PATH).drop_duplicates("compound_id").reset_index(drop=True)

    results = []
    total = len(df)
    
    for i, row in df.iterrows():
        cid = row["compound_id"]
        print(f"[{i+1}/{total}] Docking {cid}...", end=" ", flush=True)
        
        ligand_pdbqt = TEMP_DIR / f"{cid}_ligand.pdbqt"
        out_pdbqt = TEMP_DIR / f"{cid}_docked.pdbqt"
        log_file = TEMP_DIR / f"{cid}_vina.log"
        
        record = {
            "compound_id": cid,
            "docking_missing": True,
            "docking_affinity_kcal_per_mol": None,
            "docking_score": 0.0,
            "docking_pose_status": "failed",
            "failure_reason": ""
        }
        
        if not prepare_ligand(row["canonical_smiles"], ligand_pdbqt):
            record["failure_reason"] = "ligand_prep_failed"
            print("FAILED (prep)")
            results.append(record)
            continue
            
        # Construct Vina command (removed --log as it is unsupported in 1.2.x)
        cmd = [
            str(VINA_EXE),
            "--receptor", str(RECEPTOR_PATH),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(config["grid_center"][0]),
            "--center_y", str(config["grid_center"][1]),
            "--center_z", str(config["grid_center"][2]),
            "--size_x", str(config["grid_size"][0]),
            "--size_y", str(config["grid_size"][1]),
            "--size_z", str(config["grid_size"][2]),
            "--exhaustiveness", str(config["vina_exhaustiveness"]),
            "--num_modes", str(config["vina_num_modes"]),
            "--seed", str(config["vina_seed"]),
            "--out", str(out_pdbqt)
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            # Write Vina's stdout to our log file
            log_file.write_text(result.stdout)
            
            if result.returncode != 0:
                record["failure_reason"] = f"vina_error: {result.stderr[:100]}"
                print("FAILED (vina error)")
                # If it fails, save and break so you can see the error
                results.append(record)
                break
            else:
                affinity = parse_vina_affinity(out_pdbqt)
                if affinity is None:
                    record["failure_reason"] = "parse_failed"
                    print("FAILED (parse)")
                else:
                    record["docking_missing"] = False
                    record["docking_affinity_kcal_per_mol"] = affinity
                    # Locked transform: clip((-affinity - 5) / 6, 0, 1)
                    raw_score = (-affinity - 5.0) / 6.0
                    record["docking_score"] = float(np.clip(raw_score, 0.0, 1.0))
                    record["docking_pose_status"] = "top_pose_only"
                    print(f"OK (affinity: {affinity:.2f})")
                    
        except Exception as e:
            record["failure_reason"] = str(e)
            print(f"FAILED (exception)")
            
        results.append(record)

    # Save results
    out_df = pd.DataFrame(results)
    out_df.to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved {len(out_df)} docking results to {OUTPUT_PATH}")
    
    # --- NEW CODE: Save persistent poses for ChimeraX rendering ---
    POSES_DIR = ROOT / "data" / "docking_poses"
    POSES_DIR.mkdir(parents=True, exist_ok=True)
    for pdbqt_file in TEMP_DIR.glob("*_docked.pdbqt"):
        shutil.copy(pdbqt_file, POSES_DIR / pdbqt_file.name)
    print(f"Saved persistent poses to {POSES_DIR}")
    # -------------------------------------------------------------
    
    # Summary (safely handle cases where all runs failed and df is empty)
    if len(out_df) > 0 and "docking_missing" in out_df.columns:
        ok = out_df[~out_df["docking_missing"]]
        print(f"Successfully docked: {len(ok)}/{total}")
        if len(ok) > 0:
            print(f"Best affinity: {ok['docking_affinity_kcal_per_mol'].min():.2f} kcal/mol")
    else:
        print("No docking results generated.")

if __name__ == "__main__":
    main()