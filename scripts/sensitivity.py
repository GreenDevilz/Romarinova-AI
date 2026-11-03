import pandas as pd
import numpy as np
from pathlib import Path

def run_sensitivity():
    print("Running sensitivity analysis...")
    BASE_DIR = Path(__file__).resolve().parent.parent
    
    # Load the annotated screening results
    df = pd.read_csv(BASE_DIR / "results" / "screened_compounds_annotated.csv")
    
    # Apply the duplication cap (max 2 per compound)
    locked = df.sort_values("final_score", ascending=False).groupby("compound_id").head(2).head(10)
    ligand = df.sort_values("ligand_only_score", ascending=False).groupby("compound_id").head(2).head(10)
    
    # Calculate overlaps
    locked_pairs = set(zip(locked["compound_id"], locked["cultivar_source"]))
    lig_pairs = set(zip(ligand["compound_id"], ligand["cultivar_source"]))
    pair_overlap = len(locked_pairs.intersection(lig_pairs))
    
    locked_compounds = set(locked["compound_id"])
    lig_compounds = set(ligand["compound_id"])
    compound_overlap = len(locked_compounds.intersection(lig_compounds))
    
    # Determine label
    N_star = len(locked)
    U_star = len(locked_compounds)
    
    if N_star == 0:
        label = "not_applicable_NO"
    elif pair_overlap >= np.ceil(0.67 * N_star):
        label = "robust"
    elif pair_overlap >= np.ceil(0.50 * N_star):
        label = "partially_robust"
    else:
        label = "weight_sensitive"
        
    print(f"Sensitivity Label: {label} (Overlap: {pair_overlap}/{N_star})")
    
    # Write to diagnostics
    diag_path = BASE_DIR / "results" / "diagnostics.csv"
    diag_df = pd.read_csv(diag_path)
    
    # Remove old sensitivity rows if they exist
    diag_df = diag_df[diag_df["section"] != "sensitivity"]
    
    new_row = pd.DataFrame([{
        "section": "sensitivity",
        "row_id": 1,
        "pair_overlap": pair_overlap,
        "sensitivity_label": label,
        "N_star": N_star,
        "U_star": U_star,
        "notes": f"Compound overlap: {compound_overlap}"
    }])
    
    diag_df = pd.concat([diag_df, new_row], ignore_index=True)
    diag_df.to_csv(diag_path, index=False)
    print("Sensitivity results appended to diagnostics.csv")

if __name__ == "__main__":
    run_sensitivity()