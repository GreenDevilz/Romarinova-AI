"""
Post-scoring metadata fix.

Run after the scoring/shortlist notebook regenerates primary_shortlist.csv.
Idempotent. Fixes:

  1. primary_shortlist_compound_level.csv -> U rows (shortlist compounds only)
  2. shortlist_metadata.json              -> U*, screened_compound_n, all counts at
                                             unique-shortlist-compound level
  3. diagnostics.csv                      -> ledger scaffold counts refreshed from
                                             current split; row_ids renumbered;
                                             filtered_but_interesting rebuilt from
                                             screened_compounds_annotated.csv;
                                             sensitivity rows deduped to one

Cwd-independent. Run from anywhere:
    python scripts/fix_shortlist_metadata.py
"""
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

# --- Load
shortlist = pd.read_csv(RESULTS / "primary_shortlist.csv")
screened = pd.read_csv(RESULTS / "screened_compounds_annotated.csv")
compound_level_full = pd.read_csv(RESULTS / "primary_shortlist_compound_level.csv")
diag = pd.read_csv(RESULTS / "diagnostics.csv")

unique_ids = shortlist["compound_id"].unique().tolist()
pairs_per_compound = shortlist.groupby("compound_id").size().to_dict()

print(f"Shortlist: {len(shortlist)} pairs, {len(unique_ids)} unique compounds")
print(f"Screened universe: {screened['compound_id'].nunique()} compounds")

# ==================================================================
# 1. primary_shortlist_compound_level.csv -> U rows
# ==================================================================
identified = screened[screened["cultivar_evidence_level"] == "identified"]

all_cultivars_map = (
    screened.groupby("compound_id")["cultivar_source"]
    .apply(lambda s: "; ".join(sorted(s.unique())))
    .to_dict()
)
n_cultivars_map = (
    screened.groupby("compound_id")["cultivar_source"].nunique().to_dict()
)
primary_cultivar_map = (
    identified.groupby("compound_id")["cultivar_source"]
    .apply(lambda s: sorted(s.unique())[0])
    .to_dict()
)

compound_level = (
    compound_level_full[compound_level_full["compound_id"].isin(unique_ids)]
    .drop_duplicates("compound_id")
    .reset_index(drop=True)
    .copy()
)
if len(compound_level) != len(unique_ids):
    missing = set(unique_ids) - set(compound_level["compound_id"])
    print(f"WARNING: {len(compound_level)} rows vs {len(unique_ids)} expected. Missing: {missing}")

compound_level["all_cultivars"]      = compound_level["compound_id"].map(all_cultivars_map)
compound_level["n_cultivar_entries"] = compound_level["compound_id"].map(n_cultivars_map).astype(int)
compound_level["primary_cultivar"]   = compound_level["compound_id"].map(primary_cultivar_map)
compound_level["n_shortlist_pairs"]  = compound_level["compound_id"].map(pairs_per_compound).astype(int)

compound_level.to_csv(RESULTS / "primary_shortlist_compound_level.csv", index=False)
print(f"[1] compound_level: {len(compound_level)} rows")

# ==================================================================
# 2. shortlist_metadata.json
# ==================================================================
dedup = shortlist.drop_duplicates("compound_id")
metadata = {
    "N_star":                         int(len(shortlist)),
    "U_star":                         int(len(unique_ids)),
    "shortlist_pair_n":               int(len(shortlist)),
    "screened_compound_n":            int(screened["compound_id"].nunique()),
    "known_active_count":             int(dedup["known_active"].sum()),
    "high_ml_ad_extrapolation_count": int(dedup["high_ml_ad_extrapolation"].sum()),
    "ad_extrapolation_count":         int((dedup["ad_flag"] == "extrapolation").sum()),
    "ad_missing_count":               0,
    "docking_missing_count":          0,
}
(RESULTS / "shortlist_metadata.json").write_text(json.dumps(metadata, indent=2))
print(f"[2] metadata: {json.dumps(metadata)}")

# ==================================================================
# 3. diagnostics.csv
# ==================================================================
# Drop blank rows and empty-section rows
diag = diag.dropna(how="all")
diag = diag[diag["section"].notna() & (diag["section"].astype(str).str.strip() != "")]

# --- 3a. Ledger: refresh scaffold stats from current split, renumber
splits = json.loads((ROOT / "splits" / "split_indices.json").read_text())
src = pd.read_csv(ROOT / splits["curated_csv"]).reset_index(drop=True)

scaf_by_idx = {}
for i, smi in enumerate(src["canonical_smiles"]):
    m = Chem.MolFromSmiles(smi)
    scaf_by_idx[i] = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False) if m else ""

def scaffold_stats(idxs):
    counts = Counter(scaf_by_idx[i] for i in idxs)
    return len(counts), (max(counts.values()) if counts else 0)

tr_g, tr_max   = scaffold_stats(splits["indices"]["train"])
val_g, val_max = scaffold_stats(splits["indices"]["validation"])
te_g, te_max   = scaffold_stats(splits["indices"]["test"])

ledger = diag[diag["section"] == "data_feasibility_ledger"].copy()
override = {
    "scaffold_groups_train":        tr_g,
    "scaffold_groups_val":          val_g,
    "scaffold_groups_test":         te_g,
    "largest_scaffold_group_train": tr_max,
    "largest_scaffold_group_val":   val_max,
    "largest_scaffold_group_test":  te_max,
    "split_repairs_n":              len(splits.get("split_repairs", [])),
}
for metric, value in override.items():
    ledger.loc[ledger["metric"] == metric, "value"] = value

ledger = ledger.sort_values(["row_id", "metric"]).reset_index(drop=True)
ledger["row_id"] = range(1, len(ledger) + 1)

# --- 3b. filtered_but_interesting: rebuild from screened universe
not_in_shortlist = screened[~screened["compound_id"].isin(unique_ids)]
failing = (
    not_in_shortlist[not_in_shortlist["developability_score"] < 1.0]
    .drop_duplicates("compound_id")
    .sort_values("final_score", ascending=False)
    .head(5)
)
filtered_rows = []
for i, (_, row) in enumerate(failing.iterrows(), 1):
    failed = []
    if not row["lipinski"]:   failed.append("Lipinski")
    if not row["veber"]:      failed.append("Veber")
    if not row["pains"]:      failed.append("PAINS")
    if not row["aggregator"]: failed.append("aggregator")
    filtered_rows.append({
        "section":       "filtered_but_interesting",
        "row_id":        i,
        "compound_id":   row["compound_id"],
        "compound_name": row["compound_name"],
        "reason":        "Failed developability",
        "final_score":   float(row["final_score"]),
        "failed_rules":  "; ".join(failed),
    })
filtered = pd.DataFrame(filtered_rows)
print(f"[3b] filtered_but_interesting: {len(filtered)} rows")

# --- 3c. sensitivity: dedupe to one row, correct U*
sens_all = diag[diag["section"] == "sensitivity"].copy()
sens_all = sens_all[sens_all["pair_overlap"].notna()]

if len(sens_all) == 0:
    print("[3c] WARNING: no sensitivity row found; run the scoring notebook first")
    sensitivity = pd.DataFrame()
else:
    sens_all["_has_slash"] = sens_all["notes"].astype(str).str.contains(r"\d+/\d+", regex=True, na=False)
    sens_all = sens_all.sort_values("_has_slash", ascending=False)
    source = sens_all.iloc[0]
    note = str(source["notes"])
    m = re.search(r"(\d+)/\d+", note)
    compound_overlap = m.group(1) if m else re.search(r"(\d+)", note).group(1) if re.search(r"\d+", note) else "?"
    sensitivity = pd.DataFrame([{
        "section":           "sensitivity",
        "row_id":            1,
        "pair_overlap":      source["pair_overlap"],
        "sensitivity_label": source["sensitivity_label"],
        "N_star":            float(len(shortlist)),
        "U_star":            float(len(unique_ids)),
        "notes":             f"Compound overlap: {compound_overlap}/{len(unique_ids)}",
    }])
    print(f"[3c] sensitivity: {len(sens_all)} rows -> 1 "
          f"(pair_overlap={source['pair_overlap']}, compound={compound_overlap}/{len(unique_ids)})")

disagreement = diag[diag["section"] == "disagreement_patterns"].copy()

cols = [
    "section", "row_id", "metric", "split", "value", "rank", "compound_id",
    "compound_name", "cultivar_source", "source_publication", "source_attempted",
    "structure_status", "reason", "final_score", "ligand_only_score",
    "ligand_only_rank", "failed_rules", "pattern", "pair_overlap",
    "sensitivity_label", "N_star", "U_star", "notes",
]
diag_out = (
    pd.concat([filtered, disagreement, sensitivity, ledger], ignore_index=True)
    .reindex(columns=cols)
    .fillna("")
)
diag_out.to_csv(RESULTS / "diagnostics.csv", index=False)
print(f"[3] diagnostics.csv: {len(diag_out)} rows total")