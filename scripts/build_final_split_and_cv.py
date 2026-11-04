"""
Consolidated split + CV builder.
Restores the document-aware split (128/30/21, 14 train actives)
and applies class-aware CV folds on top of it.
"""
import json
from pathlib import Path
from collections import Counter

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
SP_PATH = ROOT / "splits" / "split_indices.json"
CV_PATH = ROOT / "splits" / "cv_scaffold_groups.json"

# 1. Load eligible data
src = pd.read_csv(ROOT / "data" / "curated_chembl_eligible.csv").reset_index(drop=True)
src["document_chembl_id"] = src["document_chembl_id"].fillna("__no_doc__")

# 2. Document-aware split
doc_to_idx = {}
for i, doc in enumerate(src["document_chembl_id"]):
    doc_to_idx.setdefault(doc, []).append(i)

docs = []
for doc, idxs in doc_to_idx.items():
    n_act = int((src.iloc[idxs]["activity_class"] == "active").sum())
    docs.append((doc, idxs, len(idxs), n_act))
docs.sort(key=lambda t: (-t[2], t[0])) # sort by size desc

n_total = len(src)
train_target, val_target = 0.70 * n_total, 0.85 * n_total

train_idx, val_idx, test_idx = [], [], []
for doc, idxs, size, _ in docs:
    if len(train_idx) + size <= train_target + 3:
        train_idx.extend(idxs)
    elif len(train_idx) + len(val_idx) + size <= val_target + 3:
        val_idx.extend(idxs)
    else:
        test_idx.extend(idxs)

# 3. Class-aware CV folds on training split
def scaffold_of(smi):
    m = Chem.MolFromSmiles(smi)
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False) if m else ""

scaf_by_idx = {i: scaffold_of(src.iloc[i]["canonical_smiles"]) for i in range(len(src))}
def n_act_in(idxs):
    return int((src.iloc[idxs]["activity_class"] == "active").sum())

train_groups_map = {}
for i in train_idx:
    train_groups_map.setdefault(scaf_by_idx[i], []).append(i)

active_groups = [(s, idxs) for s, idxs in train_groups_map.items() if n_act_in(idxs) > 0]
inactive_groups = [(s, idxs) for s, idxs in train_groups_map.items() if n_act_in(idxs) == 0]
active_groups.sort(key=lambda kv: (-n_act_in(kv[1]), -len(kv[1]), kv[0]))
inactive_groups.sort(key=lambda kv: (-len(kv[1]), kv[0]))

K = min(4, len(active_groups)) if active_groups else 5
folds = {str(i): [] for i in range(1, K + 1)}
fold_active_n = {str(i): 0 for i in range(1, K + 1)}
fold_compound_n = {str(i): 0 for i in range(1, K + 1)}

for scaf, idxs in active_groups:
    target = min(folds.keys(), key=lambda k: (fold_active_n[k], fold_compound_n[k], k))
    folds[target].extend(idxs)
    fold_active_n[target] += n_act_in(idxs)
    fold_compound_n[target] += len(idxs)

for scaf, idxs in inactive_groups:
    target = min(folds.keys(), key=lambda k: (fold_compound_n[k], k))
    folds[target].extend(idxs)
    fold_compound_n[target] += len(idxs)

# 4. Save
sp = {
    "schema_version": "1.2.0",
    "seed": 42,
    "curated_csv": "data/curated_chembl_eligible.csv",
    "scaffold_definition": "bemis_murcko_canonical_no_chirality",
    "split_rule": "document_aware_greedy_70_15_15",
    "n_compounds": {"train": len(train_idx), "validation": len(val_idx), "test": len(test_idx)},
    "n_actives": {"train": n_act_in(train_idx), "validation": n_act_in(val_idx), "test": n_act_in(test_idx)},
    "indices": {"train": train_idx, "validation": val_idx, "test": test_idx},
}
SP_PATH.write_text(json.dumps(sp, indent=2))

cv = {
    "schema_version": "1.2.0",
    "seed": 42,
    "n_folds": K,
    "fold_assignment_rule": "class_aware_greedy_by_active_count",
    "fold_active_n": fold_active_n,
    "fold_compound_n": fold_compound_n,
    "scaffold_of_compound": {str(i): scaf_by_idx[i] for i in train_idx},
    "folds": folds,
}
CV_PATH.write_text(json.dumps(cv, indent=2))

print(f"Split: train={len(train_idx)} ({n_act_in(train_idx)} actives), "
      f"val={len(val_idx)} ({n_act_in(val_idx)} actives), "
      f"test={len(test_idx)} ({n_act_in(test_idx)} actives)")
print(f"CV folds: {K}")
for f in sorted(folds.keys(), key=int):
    print(f"  fold {f}: {fold_compound_n[f]} compounds, {fold_active_n[f]} actives")