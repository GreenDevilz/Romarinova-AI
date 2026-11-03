"""
Document-aware re-split + diagnostic.

Satisfies protocol §5.4's intent -- no shared document_chembl_id between
train and any holdout split -- where a scaffold-first split cannot, because
a small number of large documents dominate the eligible dataset.

Usage (from repo root or anywhere):
    python scripts/apply_document_aware_split.py

The script:
  1. Prints document-size diagnostics (largest 10, largest-doc share).
  2. Builds a document-aware 70/15/15 split (whole documents per split).
  3. Rebuilds scaffold-grouped CV folds inside the training split.
  4. Verifies zero shared documents across train/val/test.

Aborts if the split cannot hit 55-80% train share (indicating a document
is too dominant for document-aware splitting).
"""
import json
import shutil
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# ---- Path resolution: repo root regardless of cwd
ROOT = Path(__file__).resolve().parents[1]
SP_PATH = ROOT / "splits" / "split_indices.json"
CV_PATH = ROOT / "splits" / "cv_scaffold_groups.json"

sp = json.loads(SP_PATH.read_text())

# Locate the eligible file, tolerate either naming (post-rename or pre-rename)
eligible_path = ROOT / sp.get("curated_csv", "data/curated_chembl_eligible.csv")
if not eligible_path.exists():
    fallback = ROOT / "data" / "curated_chembl_eligible.csv"
    if fallback.exists():
        eligible_path = fallback
    else:
        raise SystemExit(
            f"Cannot find curated eligible file. Looked at:\n"
            f"  {eligible_path}\n  {fallback}\n"
            f"Run 01_data_curation.ipynb or verify data/ contents."
        )

src = pd.read_csv(eligible_path).reset_index(drop=True)
src["document_chembl_id"] = src["document_chembl_id"].fillna("__no_doc__")
print(f"Loaded {len(src)} eligible compounds from {eligible_path.name}")

# ---------------------------------------------------------------
# 1. Diagnostic
# ---------------------------------------------------------------
doc_sizes = (
    src.groupby("document_chembl_id")
    .agg(n_compounds=("activity_class", "size"),
         n_actives=("activity_class", lambda s: (s == "active").sum()))
    .sort_values("n_compounds", ascending=False)
)

print(f"\nTotal documents: {len(doc_sizes)}")
print("Largest 10 documents:")
print(doc_sizes.head(10).to_string())

largest = doc_sizes.iloc[0]
largest_share = 100.0 * largest["n_compounds"] / len(src)
print(f"\nLargest doc: {largest.name} -- "
      f"{largest['n_compounds']} compounds ({largest_share:.1f}%), "
      f"{largest['n_actives']} actives")

# ---------------------------------------------------------------
# 2. Backup once
# ---------------------------------------------------------------
bak_sp = ROOT / "splits" / "split_indices.scaffold_backup.json"
bak_cv = ROOT / "splits" / "cv_scaffold_groups.scaffold_backup.json"
if not bak_sp.exists():
    shutil.copy(SP_PATH, bak_sp)
    shutil.copy(CV_PATH, bak_cv)
    print(f"\nBacked up scaffold-split files.")

# ---------------------------------------------------------------
# 3. Group indices by document
# ---------------------------------------------------------------
doc_to_idx = {}
for i, doc in enumerate(src["document_chembl_id"]):
    doc_to_idx.setdefault(doc, []).append(i)

docs = []
for doc, idxs in doc_to_idx.items():
    n_act = int((src.iloc[idxs]["activity_class"] == "active").sum())
    docs.append((doc, idxs, len(idxs), n_act))
docs.sort(key=lambda t: (-t[2], t[0]))

# ---------------------------------------------------------------
# 4. Greedy whole-document assignment to 70/15/15
# ---------------------------------------------------------------
n_total = len(src)
train_target = 0.70 * n_total
val_target   = 0.85 * n_total

train_idx, val_idx, test_idx = [], [], []
for doc, idxs, size, n_act in docs:
    if len(train_idx) + size <= train_target + 3:
        train_idx.extend(idxs)
    elif len(train_idx) + len(val_idx) + size <= val_target + 3:
        val_idx.extend(idxs)
    else:
        test_idx.extend(idxs)

# ---------------------------------------------------------------
# 5. Active-stratification repair (bounded, recorded)
# ---------------------------------------------------------------
def n_actives(idxs):
    return int((src.iloc[idxs]["activity_class"] == "active").sum()) if idxs else 0

repairs = []
for _ in range(3):
    if n_actives(val_idx) > 0 and n_actives(test_idx) > 0:
        break
    deficient = "test" if n_actives(test_idx) == 0 else "val"
    candidates = [(doc, idxs) for doc, idxs, size, n_act in docs
                  if set(idxs).issubset(set(train_idx)) and n_act > 0 and size <= 8]
    if not candidates:
        break
    candidates.sort(key=lambda t: (len(t[1]), t[0]))
    doc, idxs = candidates[0]
    (test_idx if deficient == "test" else val_idx).extend(idxs)
    train_idx = [i for i in train_idx if i not in set(idxs)]
    repairs.append({"document": doc, "moved_to": deficient, "size": len(idxs)})

# ---------------------------------------------------------------
# 6. Guard: split share sanity
# ---------------------------------------------------------------
pct_train = 100.0 * len(train_idx) / n_total
pct_val   = 100.0 * len(val_idx) / n_total
pct_test  = 100.0 * len(test_idx) / n_total
print(f"\nProposed split:")
print(f"  train: {len(train_idx):3d} ({pct_train:.1f}%) "
      f"-- {n_actives(train_idx)} actives")
print(f"  val:   {len(val_idx):3d} ({pct_val:.1f}%) "
      f"-- {n_actives(val_idx)} actives")
print(f"  test:  {len(test_idx):3d} ({pct_test:.1f}%) "
      f"-- {n_actives(test_idx)} actives")
print(f"  repairs: {len(repairs)}")

if not (55 <= pct_train <= 80):
    raise SystemExit(
        f"ABORT: train share {pct_train:.1f}% outside 55-80% band. "
        f"Largest document is too dominant for document-aware splitting "
        f"({largest_share:.1f}% of eligible). Switch to Path B: keep the "
        f"scaffold split, report leakage as a documented caveat."
    )

# ---------------------------------------------------------------
# 7. Rebuild scaffold-grouped CV folds within training split
# ---------------------------------------------------------------
def scaffold_of(smi):
    m = Chem.MolFromSmiles(smi)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False) if m else ""

scaf_by_idx = {i: scaffold_of(src.iloc[i]["canonical_smiles"]) for i in train_idx}
groups = {}
for i in train_idx:
    groups.setdefault(scaf_by_idx[i], []).append(i)
groups_sorted = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))

folds = {str(i): [] for i in range(1, 6)}
for i, (_, idxs) in enumerate(groups_sorted):
    folds[str((i % 5) + 1)].extend(idxs)

# ---------------------------------------------------------------
# 8. Persist
# ---------------------------------------------------------------
sp["indices"] = {
    "train":      train_idx,
    "validation": val_idx,
    "test":       test_idx,
}
sp["n_compounds"] = {"train": len(train_idx), "validation": len(val_idx), "test": len(test_idx)}
sp["n_actives"]   = {
    "train":      n_actives(train_idx),
    "validation": n_actives(val_idx),
    "test":       n_actives(test_idx),
}
sp["split_rule"] = "document_aware_greedy_70_15_15_with_active_repair"
sp["split_repairs"] = repairs
sp["document_aware"] = True
sp["document_aware_rationale"] = (
    "Scaffold-first split combined with document exclusion reduced training "
    "to 30/179 (16.8%). Documents are now treated as atomic assignment units; "
    "scaffold-grouped CV within training preserves scaffold-aware model selection."
)
SP_PATH.write_text(json.dumps(sp, indent=2))

cv = json.loads(CV_PATH.read_text())
cv["folds"] = folds
cv["fold_assignment_rule"] = "size_desc_round_robin_over_train_scaffold_groups"
CV_PATH.write_text(json.dumps(cv, indent=2))

# ---------------------------------------------------------------
# 9. Verification
# ---------------------------------------------------------------
def documents(idxs):
    return set(src.iloc[idxs]["document_chembl_id"].dropna())

td, vd, sd = documents(train_idx), documents(val_idx), documents(test_idx)
print(f"\n=== VERIFICATION ===")
print(f"train ∩ val  docs: {len(td & vd)}  (target 0)")
print(f"train ∩ test docs: {len(td & sd)}  (target 0)")
print(f"val   ∩ test docs: {len(vd & sd)}  (target 0)")
print(f"\nOK -- split persisted to {SP_PATH.name}")
print(f"Next: rerun notebooks/02_stage1_training.ipynb end-to-end.")