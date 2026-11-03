"""
Validate web/data_contract.json against itself and against the source CSVs.

Checks:
  * Top-level and metadata keys present.
  * Row counts match metadata (N_star, U_star, screened_compound_n).
  * shortlist_pairs and compound_level_view reference the same compounds.
  * compound_level_view ⊆ screened_compounds (referential integrity).
  * Flag counts consistent (ad_extrapolation_count, high_ml_ad_extrapolation_count).

Exit 0 on pass, 1 on failure.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "web" / "data_contract.json"


def main():
    print(f"Validating {CONTRACT}...")

    if not CONTRACT.exists():
        print("FAIL: data_contract.json does not exist.")
        return 1
    try:
        data = json.loads(CONTRACT.read_text())
    except json.JSONDecodeError as e:
        print(f"FAIL: not valid JSON: {e}")
        return 1

    errors = []

    for key in ["metadata", "shortlist_pairs", "compound_level_view",
                "screened_compounds", "diagnostics"]:
        if key not in data:
            errors.append(f"missing top-level key '{key}'")

    meta = data.get("metadata", {})
    for key in ["N_star", "U_star", "screened_compound_n",
                "known_active_count", "high_ml_ad_extrapolation_count",
                "ad_extrapolation_count", "ad_missing_count"]:
        if key not in meta:
            errors.append(f"missing metadata key '{key}'")

    if errors:
        _report(errors)
        return 1

    sl  = data["shortlist_pairs"]
    clv = data["compound_level_view"]
    sc  = data["screened_compounds"]

    # --- Counts
    if len(sl) != meta["N_star"]:
        errors.append(f"shortlist_pairs has {len(sl)} rows; N_star = {meta['N_star']}")
    if len(clv) != meta["U_star"]:
        errors.append(f"compound_level_view has {len(clv)} rows; U_star = {meta['U_star']}")
    if len(sc) != meta["screened_compound_n"]:
        errors.append(f"screened_compounds has {len(sc)} rows; "
                      f"screened_compound_n = {meta['screened_compound_n']}")

    # --- Referential integrity
    sl_ids  = {r["compound_id"] for r in sl}
    clv_ids = {r["compound_id"] for r in clv}
    sc_ids  = {r["compound_id"] for r in sc}

    if clv_ids - sl_ids:
        errors.append(f"compound_level_view has compounds absent from shortlist: "
                      f"{sorted(clv_ids - sl_ids)}")
    if sl_ids - clv_ids:
        errors.append(f"shortlist has compounds absent from compound_level_view: "
                      f"{sorted(sl_ids - clv_ids)}")
    if clv_ids - sc_ids:
        errors.append(f"compound_level_view has compounds absent from screened_compounds: "
                      f"{sorted(clv_ids - sc_ids)}")

    # --- Compound-level uniqueness
    if len(clv_ids) != len(clv):
        errors.append(f"compound_level_view has "
                      f"{len(clv) - len(clv_ids)} duplicate compound_id(s)")

    # --- Flag consistency
    extrapolation_n = sum(1 for r in clv if r.get("ad_flag") == "extrapolation")
    if extrapolation_n != meta["ad_extrapolation_count"]:
        errors.append(f"compound_level_view has {extrapolation_n} ad_flag=='extrapolation'; "
                      f"metadata ad_extrapolation_count = {meta['ad_extrapolation_count']}")

    high_ml_ad_n = sum(1 for r in clv if r.get("high_ml_ad_extrapolation") is True)
    if high_ml_ad_n != meta["high_ml_ad_extrapolation_count"]:
        errors.append(f"compound_level_view has {high_ml_ad_n} high_ml_ad_extrapolation==True; "
                      f"metadata high_ml_ad_extrapolation_count = {meta['high_ml_ad_extrapolation_count']}")

    if errors:
        _report(errors)
        return 1

    print(f"OK: {meta['N_star']} pairs / {meta['U_star']} compounds / "
          f"{meta['screened_compound_n']} screened / "
          f"{len(data['diagnostics'])} diagnostics rows")
    return 0


def _report(errors):
    print("FAIL")
    for e in errors:
        print(f"  X {e}")


if __name__ == "__main__":
    sys.exit(main())