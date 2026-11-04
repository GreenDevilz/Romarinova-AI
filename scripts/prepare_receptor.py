import os
import urllib.request
import subprocess
import json
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PDB_ID = "3L4W"
PDB_PATH = ROOT / "data" / f"{PDB_ID}.pdb"
PDBQT_PATH = ROOT / "data" / f"{PDB_ID}_receptor.pdbqt"
CONFIG_PATH = ROOT / "models" / "docking_config.json"

def main():
    if not PDB_PATH.exists():
        print(f"Downloading {PDB_ID} from RCSB...")
        url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
        urllib.request.urlretrieve(url, PDB_PATH)
        print("Download complete.")

    # 1. Extract ligand centroid to set the grid box exactly over the active site
    coords = []
    with open(PDB_PATH, "r") as f:
        for line in f:
            if line.startswith("HETATM") and "HOH" not in line:
                try:
                    # Extract x, y, z from PDB fixed-width columns
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x,y,z])
                except ValueError:
                    continue

    if not coords:
        raise ValueError("No non-water HETATM records found. Cannot set grid center.")
    
    coords = np.array(coords)
    center = coords.mean(axis=0)
    print(f"Calculated true ligand center: {center}")

    # 2. Update docking_config.json with the calculated center
    config = json.loads(CONFIG_PATH.read_text())
    config["grid_center"] = [round(c, 2) for c in center]
    CONFIG_PATH.write_text(json.dumps(config, indent=2))
    print(f"Updated {CONFIG_PATH} with new grid_center.")

    # 3. Prepare receptor with OpenBabel (strip waters first)
    print("Preparing receptor with OpenBabel...")
    clean_pdb = ROOT / "data" / f"{PDB_ID}_clean.pdb"
    with open(PDB_PATH, "r") as fin, open(clean_pdb, "w") as fout:
        for line in fin:
            # Remove crystallographic waters
            if line.startswith("HETATM") and "HOH" in line:
                continue
            fout.write(line)
            
    cmd = [
        "obabel", 
        str(clean_pdb), 
        "-O", str(PDBQT_PATH),
        "-xr", "-p", "7.4"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"OpenBabel failed:\n{result.stderr}")
    
    # Clean up temporary PDB
    clean_pdb.unlink()
    
    print(f"Receptor saved to {PDBQT_PATH}")

if __name__ == "__main__":
    main()