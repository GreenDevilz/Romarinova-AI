"""
Post-rerun verification.

Run after:
  * 02_stage1_training.ipynb           (retrain on document-aware split)
  * 04_stage2_similarity.ipynb         (fill similarity-only row)
  * 05_stage3_developability.ipynb     (real rules.json)
  * scoring / shortlist notebook       (regenerate shortlist)
  * Batch 2 metadata fix cell          (fix U*, screened_compound_n, diagnostics)

Usage from anywhere:
    python scripts/verify_post_rerun.py

This script does NOT modify any file. It reads artifacts and prints a
consolidated report to copy into the markdown reconciliation pass.
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

print("=" * 60)
print("POST-RERUN VERIFICATION")
print("=" * 60)

# ---------------------------------------------------------------
# CHECK 1 -- Applicability domain status of the 28 rosemary compounds
# ---------------------------------------------------------------
print("\n[CHECK 1] Applicability Domain")
ad_path = RESULTS / "ad_results.csv"
if not ad_path.exists():
    print(f"  MISSING: {ad_path}")
else:
    ad = pd.read_csv(ad_path)
    counts = ad["ad_flag"].value_counts().to_dict()
    print(f"  Total rosemary compounds: {len(ad)}")
    print(f"  Flag distribution: {counts}")
    print(f"  ad_distance range: {ad['ad_distance'].min():.4f} – {ad['ad_distance'].max():.4f}")

    # Threshold: check model_config first, then fall back
    mc_path = ROOT / "models" / "model_config.json"
    if mc_path.exists():
        mc = json.loads(mc_path.read_text())
        thr = mc.get("ad_threshold") or mc.get("ad_loo_threshold")
        if thr:
            print(f"  AD threshold (model_config.json): {thr:.4f}")
        else:
            print(f"  AD threshold: not stored in model_config.json")
    else:
        print(f"  model_config.json not found")

# ---------------------------------------------------------------
# CHECK 2 -- Shortlist composition and U*
# ---------------------------------------------------------------
print("\n[CHECK 2] Shortlist composition")
sl_path = RESULTS / "primary_shortlist.csv"
if not sl_path.exists():
    print(f"  MISSING: {sl_path}")
else:
    sl = pd.read_csv(sl_path)
    print(f"  Pairs: {len(sl)}")
    print(f"  Unique compounds: {sl['compound_id'].nunique()}")
    print(f"\n  Compound-level summary:")
    summary = (
        sl.groupby("compound_id")
        .agg(
            n_pairs=("compound_id", "size"),
            final_rank=("final_rank", "first"),
            final_score=("final_score", "first"),
            ml_rank=("ML_rank_score", "first"),
            ad_flag=("ad_flag", "first"),
            developability=("developability_score", "first"),
        )
        .sort_values("final_rank")
    )
    print(summary.to_string())

    # Check whether characteristic diterpenes made the shortlist
    diterpenes = ["carnosic_acid", "carnosol", "rosmarinic_acid"]
    in_shortlist = [d for d in diterpenes if d in sl["compound_id"].values]
    print(f"\n  Characteristic diterpenes in shortlist: "
          f"{in_shortlist if in_shortlist else 'none (all in diagnostics)'}")

# ---------------------------------------------------------------
# CHECK 3 -- Sensitivity overlap
# ---------------------------------------------------------------
print("\n[CHECK 3] Sensitivity overlap")
diag_path = RESULTS / "diagnostics.csv"
if not diag_path.exists():
    print(f"  MISSING: {diag_path}")
else:
    diag = pd.read_csv(diag_path)
    sens = diag[diag["section"] == "sensitivity"] if "section" in diag.columns else pd.DataFrame()
    if sens.empty:
        print(f"  No sensitivity row found in diagnostics.csv")
    else:
        cols = [c for c in ["pair_overlap", "sensitivity_label", "N_star", "U_star", "notes"]
                if c in sens.columns]
        print(sens[cols].to_string(index=False))

# ---------------------------------------------------------------
# Sanity: shortlist_metadata in corrected format
# ---------------------------------------------------------------
print("\n[SANITY] shortlist_metadata.json")
sm_path = RESULTS / "shortlist_metadata.json"
if not sm_path.exists():
    print(f"  MISSING: {sm_path}")
else:
    sm = json.loads(sm_path.read_text())
    expected_keys = ["N_star", "U_star", "screened_compound_n"]
    missing = [k for k in expected_keys if k not in sm]
    if missing:
        print(f"  MISSING keys: {missing}")
        print(f"  --> Batch 2 metadata fix cell has NOT run")
    else:
        print(f"  All expected keys present")
        for k in ["N_star", "U_star", "screened_compound_n",
                  "known_active_count", "high_ml_ad_extrapolation_count",
                  "ad_extrapolation_count", "ad_missing_count"]:
            if k in sm:
                print(f"  {k}: {sm[k]}")

print("\n" + "=" * 60)
print("END OF REPORT")
print("=" * 60)