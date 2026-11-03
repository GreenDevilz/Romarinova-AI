import json
from pathlib import Path
import pandas as pd

et   = pd.read_csv("results/evaluation_table.csv")
sm   = json.loads(Path("results/shortlist_metadata.json").read_text())
sp   = json.loads(Path("splits/split_indices.json").read_text())
ad   = pd.read_csv("results/ad_results.csv")
diag = pd.read_csv("results/diagnostics.csv")

test = et[(et.model == "Selected RF") & (et.split == "test")].iloc[0]
val  = et[(et.model == "Selected RF") & (et.split == "validation")].iloc[0]
sens = diag[diag.section == "sensitivity"].iloc[0] if (diag.section == "sensitivity").any() else {}

print("=== REPORT / DEVPOST VALUES ===")
print(f"test PR-AUC      : {test.pr_auc:.3f}")
print(f"test ROC-AUC     : {test.roc_auc:.3f}")
print(f"test N_actives   : {test.test_n_actives}")
print(f"test N_inactives : {test.test_n_inactives}")
print(f"test flags       : ef_unstable={test.ef_unstable}, "
      f"ef_not_reported={test.ef_not_reported}, ci_not_reported={test.ci_not_reported}")
print(f"val PR-AUC       : {val.pr_auc:.3f}")
print(f"val ROC-AUC      : {val.roc_auc:.3f}")
print(f"train actives    : {sp['n_actives']['train']}  (of {sp['n_compounds']['train']})")
print(f"val actives      : {sp['n_actives']['validation']}  (of {sp['n_compounds']['validation']})")
print(f"test actives     : {sp['n_actives']['test']}  (of {sp['n_compounds']['test']})")
print(f"shortlist N*     : {sm['N_star']}")
print(f"shortlist U*     : {sm['U_star']}")
print(f"screened universe: {sm['screened_compound_n']}")
print(f"high_ml_ad_count : {sm['high_ml_ad_extrapolation_count']}")
print(f"ad_extrapolation : {(ad.ad_flag == 'extrapolation').sum()}/{len(ad)} rosemary")
print(f"sensitivity      : pair_overlap={sens.get('pair_overlap', '?')}, "
      f"label={sens.get('sensitivity_label', '?')}, "
      f"notes={sens.get('notes', '?')}")