import pandas as pd
import json
from pathlib import Path

print("Starting artifact generation...")

# Load the locked data
shortlist = pd.read_csv("results/primary_shortlist.csv")
compound_level = pd.read_csv("results/primary_shortlist_compound_level.csv")
annotated = pd.read_csv("results/screened_compounds_annotated.csv")
meta = json.loads(Path("results/shortlist_metadata.json").read_text())
diag = pd.read_csv("results/diagnostics.csv")

# Fill NaNs with None for clean JSON
shortlist = shortlist.where(pd.notnull(shortlist), None)
compound_level = compound_level.where(pd.notnull(compound_level), None)
annotated = annotated.where(pd.notnull(annotated), None)
diag = diag.where(pd.notnull(diag), None)

# Build the payload
data_contract = {
    "metadata": meta,
    "shortlist_pairs": shortlist.to_dict(orient="records"),
    "compound_level_view": compound_level.to_dict(orient="records"),
    "screened_compounds": annotated.to_dict(orient="records"),
    "diagnostics": {
        "ledger": diag[diag["section"] == "data_feasibility_ledger"].to_dict(orient="records"),
        "filtered_but_interesting": diag[diag["section"] == "filtered_but_interesting"].to_dict(orient="records"),
        "disagreement_patterns": diag[diag["section"] == "disagreement_patterns"].to_dict(orient="records"),
        "sensitivity": diag[diag["section"] == "sensitivity"].to_dict(orient="records")
    }
}

# Write to web folder for P3
Path("web").mkdir(exist_ok=True)
Path("web/data_contract.json").write_text(json.dumps(data_contract, indent=2))
print("Successfully wrote web/data_contract.json")

# --- Generate Figures ---
import plotly.express as px
Path("figures").mkdir(exist_ok=True)

# Figure 1: Score distribution (Locked vs Ligand-Only)
fig1 = px.scatter(
    annotated, x="ligand_only_score", y="final_score", 
    color="high_ml_ad_extrapolation", hover_name="compound_name",
    title="Locked Score vs Ligand-Only Score (Colored by AD Extrapolation)",
    labels={"ligand_only_score": "Ligand-Only Score", "final_score": "Locked Final Score"}
)
fig1.write_image("figures/sensitivity_overlap.png")
print("Wrote figures/sensitivity_overlap.png")

print("Artifact generation complete.")