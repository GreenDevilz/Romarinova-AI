import json
from pathlib import Path
import pandas as pd

print("=== CANONICAL FILES ===")
cur = pd.read_csv("data/curated_chembl.csv")
print(f"data/curated_chembl.csv: {len(cur)} rows  {cur['activity_class'].value_counts().to_dict()}")

elig = Path("data/curated_chembl_eligible.csv")
if elig.exists():
    el = pd.read_csv(elig)
    print(f"data/curated_chembl_eligible.csv: {len(el)} rows  {el['activity_class'].value_counts().to_dict()}")
else:
    print("data/curated_chembl_eligible.csv: MISSING  <-- Fix A not applied")

sp = json.loads(Path("splits/split_indices.json").read_text())
print(f"\nsplit_indices.json curated_csv field = {sp.get('curated_csv')!r}")
src = pd.read_csv(sp.get("curated_csv", "data/curated_chembl.csv"))
for name in ["train", "validation", "test"]:
    idxs = sp["indices"][name]
    cls = src.iloc[idxs]["activity_class"].value_counts().to_dict()
    print(f"  {name:12s}: {len(idxs):3d} compounds  {cls}")

print("\n=== DOCUMENT LEAKAGE ===")
def docs_of(idxs):
    return set(src.iloc[idxs]["document_chembl_id"].dropna())
td, vd, sd = docs_of(sp["indices"]["train"]), docs_of(sp["indices"]["validation"]), docs_of(sp["indices"]["test"])
print(f"train ∩ test docs: {len(td & sd)}  (protocol: 0)")
print(f"train ∩ val  docs: {len(td & vd)}  (protocol: 0)")
print(f"val   ∩ test docs: {len(vd & sd)}  (protocol does not require 0 here, but flag)")