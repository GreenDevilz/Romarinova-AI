import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

# Keys the frontend (app.py) reads from shortlist_metadata.json. If any of
# these is missing after this script runs, we raise — that's a pipeline bug,
# not a UI bug, and it must be caught before deployment.
REQUIRED_FOR_UI = {
    "N_star",
    "U_star",
    "shortlist_pair_n",
    "screened_compound_n",
    "known_active_count",
    "high_ml_ad_extrapolation_count",
    "ad_extrapolation_count",
    "ad_missing_count",
    "docking_missing_count",
    "ml_rank_top_band_value",
    "ml_rank_top_band_count",
    "ml_rank_top_band_share",
    "ml_rank_unique_values_n",
    "pair_overlap",
    "sensitivity_label",
    "coverage_gap",
}


# ==================================================================
# 0. Load
# ==================================================================
shortlist = pd.read_csv(RESULTS / "primary_shortlist.csv")
screened = pd.read_csv(RESULTS / "screened_compounds_annotated.csv")
compound_level_full = pd.read_csv(RESULTS / "primary_shortlist_compound_level.csv")
diag = pd.read_csv(RESULTS / "diagnostics.csv")

meta_path = RESULTS / "shortlist_metadata.json"
prior_meta: dict = {}
if meta_path.exists():
    try:
        prior_meta = json.loads(meta_path.read_text())
    except Exception as exc:
        print(f"WARNING: prior shortlist_metadata.json unreadable ({exc}); starting empty")
        prior_meta = {}

unique_ids = shortlist["compound_id"].unique().tolist()
pairs_per_compound = shortlist.groupby("compound_id").size().to_dict()
screened_n = int(screened["compound_id"].nunique())

print(f"Shortlist: {len(shortlist)} pairs, {len(unique_ids)} unique compounds")
print(f"Screened universe: {screened_n} compounds")
print(f"Prior metadata keys present: {sorted(prior_meta.keys())}")


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

# ---- 2a. ML tie-band stats: preserve, then re-derive anything missing ----
# Three sources of truth, in priority order:
#   (i)   prior_meta (what notebook 06 wrote last run)
#   (ii)  results/ml_ad_scores.csv (the scored universe)
#   (iii) screened_compounds_annotated.csv (last resort)
# Any key absent from (i) is recovered from (ii) or (iii).
ml_band: dict = {}
for _k in (
    "ml_rank_top_band_value",
    "ml_rank_top_band_count",
    "ml_rank_unique_values_n",
    "ml_rank_top_band_share",
):
    if prior_meta.get(_k) is not None:
        ml_band[_k] = prior_meta[_k]

_need_ml = {
    "ml_rank_top_band_value",
    "ml_rank_top_band_count",
    "ml_rank_unique_values_n",
} - ml_band.keys()

if _need_ml:
    _src = None
    _ml_path = RESULTS / "ml_ad_scores.csv"
    if _ml_path.exists():
        _src = pd.read_csv(_ml_path)
    elif "ML_rank_score" in screened.columns:
        _src = screened

    if _src is not None and "ML_rank_score" in _src.columns:
        _vals = pd.to_numeric(_src["ML_rank_score"], errors="coerce").dropna()
        if not _vals.empty:
            # Round to 6 dp before value_counts so floating-point noise in
            # otherwise-identical scores doesn't fracture the tie-band.
            _counts = _vals.round(6).value_counts()
            ml_band.setdefault("ml_rank_top_band_value", float(_counts.index[0]))
            ml_band.setdefault("ml_rank_top_band_count", int(_counts.iloc[0]))
            ml_band.setdefault("ml_rank_unique_values_n", int(len(_counts)))

# Share is always recomputed from count / screened_n — cheap, exact, and
# keeps the three ML fields internally consistent regardless of source.
if "ml_rank_top_band_count" in ml_band and screened_n:
    ml_band["ml_rank_top_band_share"] = round(
        ml_band["ml_rank_top_band_count"] / screened_n, 4
    )

print(f"[2a] ml_band: {ml_band}")

# ---- 2b. pair_overlap / sensitivity_label from diagnostics (authoritative) ----
# The scoring notebook (cell 4) writes a single row in the 'sensitivity'
# section with pair_overlap as an INTEGER count and sensitivity_label as a
# string. Reform pair_overlap to "N/N_star" for the UI here.
pair_overlap_str = "0/0"
sensitivity_label = "unknown"
if "section" in diag.columns:
    _sens = diag[diag["section"] == "sensitivity"]
    if not _sens.empty:
        _row = _sens.iloc[0]
        _po = _row.get("pair_overlap")
        _n_star = int(len(shortlist))
        if pd.notna(_po):
            try:
                pair_overlap_str = f"{int(_po)}/{_n_star}"
            except (TypeError, ValueError):
                pass
        _sl = _row.get("sensitivity_label")
        if pd.notna(_sl) and str(_sl).strip():
            sensitivity_label = str(_sl).strip()
    else:
        print("[2b] WARNING: no sensitivity row in diagnostics.csv — "
              "run the scoring notebook first")

# ---- 2c. coverage_gap from missing counts ----
ad_missing     = int(prior_meta.get("ad_missing_count", 0) or 0)
docking_missing = int(prior_meta.get("docking_missing_count", 0) or 0)
if ad_missing == 0 and docking_missing == 0:
    coverage_gap = "None"
else:
    _parts = []
    if ad_missing:      _parts.append(f"{ad_missing} AD")
    if docking_missing: _parts.append(f"{docking_missing} docking")
    coverage_gap = ", ".join(_parts) + " missing"

# ---- 2d. Assemble (explicit computed keys override stale prior values) ----
metadata = {
    "N_star":                         int(len(shortlist)),
    "U_star":                         int(len(unique_ids)),
    "shortlist_pair_n":               int(len(shortlist)),
    "screened_compound_n":            screened_n,
    "known_active_count":             int(dedup["known_active"].sum()),
    "high_ml_ad_extrapolation_count": int(dedup["high_ml_ad_extrapolation"].sum()),
    "ad_extrapolation_count":         int((dedup["ad_flag"] == "extrapolation").sum()),
    "ad_missing_count":               ad_missing,
    "docking_missing_count":          docking_missing,
    "pair_overlap":                   pair_overlap_str,
    "sensitivity_label":              sensitivity_label,
    "coverage_gap":                   coverage_gap,
}

# ML fields override anything in `metadata` (none today, but keeps intent clear).
metadata.update(ml_band)

# Preserve any extra keys the notebook wrote that this script doesn't manage
# (forward-compatibility for future fields). setdefault ensures we never
# overwrite what we computed above.
for _k, _v in prior_meta.items():
    metadata.setdefault(_k, _v)

meta_path.write_text(json.dumps(metadata, indent=2))
print(f"[2] metadata written ({len(metadata)} keys)")
for _k in sorted(metadata):
    print(f"      {_k}: {metadata[_k]!r}")


# ==================================================================
# 3. diagnostics.csv
# ==================================================================
# Drop blank rows and empty-section rows
diag = diag.dropna(how="all")
diag = diag[diag["section"].notna() & (diag["section"].astype(str).str.strip() != "")]

# ---- 3a. Ledger: refresh scaffold stats from current split, renumber ----
splits = json.loads((ROOT / "splits" / "split_indices.json").read_text())
src = pd.read_csv(ROOT / splits["curated_csv"]).reset_index(drop=True)

scaf_by_idx = {}
for i, smi in enumerate(src["canonical_smiles"]):
    m = Chem.MolFromSmiles(smi)
    scaf_by_idx[i] = (
        MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
        if m else ""
    )

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

# ---- 3b. filtered_but_interesting: rebuild from screened universe ----
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

# ---- 3c. sensitivity: dedupe to one row, correct U_star ----
sens_all = diag[diag["section"] == "sensitivity"].copy()
sens_all = sens_all[sens_all["pair_overlap"].notna()]

if len(sens_all) == 0:
    print("[3c] WARNING: no sensitivity row found; run the scoring notebook first")
    sensitivity = pd.DataFrame()
else:
    sens_all["_has_slash"] = (
        sens_all["notes"].astype(str).str.contains(r"\d+/\d+", regex=True, na=False)
    )
    sens_all = sens_all.sort_values("_has_slash", ascending=False)
    source = sens_all.iloc[0]
    note = str(source["notes"])
    m = re.search(r"(\d+)/\d+", note)
    compound_overlap = (
        m.group(1) if m
        else re.search(r"(\d+)", note).group(1) if re.search(r"\d+", note)
        else "?"
    )
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
          f"(pair_overlap={source['pair_overlap']}, "
          f"compound={compound_overlap}/{len(unique_ids)})")

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


# ==================================================================
# 4. Contract assertion — fail loudly if the UI would render '—'
# ==================================================================
_missing_from_meta = REQUIRED_FOR_UI - metadata.keys()
if _missing_from_meta:
    raise RuntimeError(
        "shortlist_metadata.json is missing UI-required keys: "
        f"{sorted(_missing_from_meta)}. "
        "Add a derivation path in fix_shortlist_metadata.py §2 before deploying."
    )

# Sanity: pair_overlap and N_star must agree on the denominator.
_n_num, _, _n_den = pair_overlap_str.partition("/")
if _n_den and int(_n_den) != int(metadata["N_star"]):
    raise RuntimeError(
        f"pair_overlap ({pair_overlap_str}) disagrees with N_star "
        f"({metadata['N_star']}) — diagnostics are stale. "
        "Rerun notebooks/06_sensitivity_analysis.ipynb first."
    )

print(f"[4] contract OK — all {len(REQUIRED_FOR_UI)} UI keys present.")