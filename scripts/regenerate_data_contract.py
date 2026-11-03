"""
Regenerate web/data_contract.json from the current results.

Fixes from previous version:
  * screened_compounds deduped to one row per compound, with a nested
    `cultivars` array carrying the pair-level provenance. Length now
    matches metadata.screened_compound_n (28), not the raw CSV length (64).
  * screened_pair_n added to metadata so both counts are explicit.
  * Referential integrity: every compound in compound_level_view is
    guaranteed to appear in screened_compounds.

Run after fix_shortlist_metadata.py:
    python scripts/regenerate_data_contract.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
WEB = ROOT / "web"
WEB.mkdir(exist_ok=True)

shortlist      = pd.read_csv(RESULTS / "primary_shortlist.csv")
compound_level = pd.read_csv(RESULTS / "primary_shortlist_compound_level.csv")
screened_raw   = pd.read_csv(RESULTS / "screened_compounds_annotated.csv")
diagnostics    = pd.read_csv(RESULTS / "diagnostics.csv")
metadata       = json.loads((RESULTS / "shortlist_metadata.json").read_text())

# --- Dedupe screened compounds to one row per compound, with nested cultivars
PAIR_LEVEL_FIELDS = [
    "cultivar_source", "cultivar_evidence_level",
    "source_publication", "source_locator", "retrieval_rank",
]
COMPOUND_LEVEL_FIELDS = [
    c for c in screened_raw.columns if c not in PAIR_LEVEL_FIELDS
]

def compound_record(compound_id, group):
    first = group.iloc[0]
    compound_row = {k: first[k] for k in COMPOUND_LEVEL_FIELDS}

    # Cast numpy types to plain Python for JSON
    for k, v in list(compound_row.items()):
        if pd.isna(v):
            compound_row[k] = None
        elif hasattr(v, "item"):
            compound_row[k] = v.item()

    compounds_cultivars = []
    for _, row in group.iterrows():
        cultivars_entry = {k: row[k] for k in PAIR_LEVEL_FIELDS}
        for k, v in list(cultivars_entry.items()):
            if pd.isna(v):
                cultivars_entry[k] = None
            elif hasattr(v, "item"):
                cultivars_entry[k] = v.item()
        compounds_cultivars.append(cultivars_entry)

    compound_row["cultivars"] = compounds_cultivars
    compound_row["n_cultivars"] = len(compounds_cultivars)
    return compound_row

screened_compounds_list = [
    compound_record(cid, group)
    for cid, group in screened_raw.groupby("compound_id", sort=False)
]

# --- Augment metadata with both counts
metadata_out = dict(metadata)
metadata_out["screened_pair_n"] = int(len(screened_raw))

# --- Build contract
def records(df):
    return json.loads(df.to_json(orient="records"))

contract = {
    "schema_version": "1.2.0",
    "generated_from": {
        "shortlist_pairs_n":          len(shortlist),
        "compound_level_view_n":      len(compound_level),
        "screened_compounds_n":       len(screened_compounds_list),
        "screened_pairs_raw_n":       len(screened_raw),
        "diagnostics_rows_n":         len(diagnostics),
    },
    "metadata":            metadata_out,
    "shortlist_pairs":     records(shortlist),
    "compound_level_view": records(compound_level),
    "screened_compounds":  screened_compounds_list,
    "diagnostics":         records(diagnostics),
}

(WEB / "data_contract.json").write_text(json.dumps(contract, indent=2))

# --- Cross-check
cl_ids = {r["compound_id"] for r in contract["compound_level_view"]}
sc_ids = {r["compound_id"] for r in contract["screened_compounds"]}
missing_in_screened = cl_ids - sc_ids
if missing_in_screened:
    raise SystemExit(
        f"FAIL: compound_level_view contains compounds absent from "
        f"screened_compounds: {sorted(missing_in_screened)}"
    )

print("web/data_contract.json regenerated:")
print(f"  shortlist_pairs     : {len(contract['shortlist_pairs'])}")
print(f"  compound_level_view : {len(contract['compound_level_view'])}")
print(f"  screened_compounds  : {len(contract['screened_compounds'])} (unique compounds)")
print(f"  screened_pairs_raw  : {len(screened_raw)} (pair-level, preserved in CSV)")
print(f"  diagnostics rows    : {len(contract['diagnostics'])}")
print(f"  metadata U_star     : {metadata_out['U_star']}")
print(f"  metadata N_star     : {metadata_out['N_star']}")
print(f"  metadata screened_n : {metadata_out['screened_compound_n']}")
print(f"  metadata screened_p : {metadata_out['screened_pair_n']}")