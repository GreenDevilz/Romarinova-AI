"""
Post-lock correction script.

Applies the protocol-compliant corrections recorded in compliance.md §4
decision-log rows 25-31:

  * Excludes training compounds whose document_chembl_id is shared with
    the validation or test split (protocol §5.4).
  * Rebuilds CV folds on the surviving training set.
  * Renames data/curated_chembl.csv to hold the eligible-only set and
    archives the full set under data/curated_chembl_with_ambiguous.csv.

Run from the repo root:
    python scripts/apply_post_lock_fixes.py

Idempotent: safe to rerun. Backs up split files on first run only.
"""
import json
import shutil
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
SP_PATH  = ROOT / "splits" / "split_indices.json"
CV_PATH  = ROOT / "splits" / "cv_scaffold_groups.json"
SP_BAK   = ROOT / "splits" / "split_indices.backup.json"
CV_BAK   = ROOT / "splits" / "cv_scaffold_groups.backup.json"
DATA_DIR = ROOT / "data"

# ---------------------------------------------------------------
# 0. Backup (first run only)
# ---------------------------------------------------------------
if not SP_BAK.exists():
    shutil.copy(SP_PATH, SP_BAK)
    shutil.copy(CV_PATH, CV_BAK)
    print(f"Backed up split files to {SP_BAK.name} / {CV_BAK.name}")

# ---------------------------------------------------------------
# 1. Document-leakage exclusion from training
# ---------------------------------------------------------------
sp  = json.loads(SP_PATH.read_text())
cv  = json.loads(CV_PATH.read_text())
src = pd.read_csv(ROOT / sp["curated_csv"]).reset_index(drop=True)

if sp.get("document_leakage_removed_from_train") is None:
    protected_docs = (
        set(src.iloc[sp["indices"]["validation"]]["document_chembl_id"].dropna())
        | set(src.iloc[sp["indices"]["test"]]["document_chembl_id"].dropna())
    )
    train_kept, train_dropped = [], []
    for i in sp["indices"]["train"]:
        if src.iloc[i]["document_chembl_id"] in protected_docs:
            train_dropped.append(int(i))
        else:
            train_kept.append(int(i))

    n_before = len(sp["indices"]["train"])
    n_after  = len(train_kept)
    pct      = 100.0 * n_after / len(src)
    print(f"Train: {n_before} -> {n_after}  "
          f"(dropped {len(train_dropped)}, {pct:.1f}% of eligible)")

    if pct < 40.0:
        raise SystemExit(
            f"ABORT: post-hoc exclusion would shrink training to {n_after} "
            f"compounds ({pct:.1f}% of eligible). Escalate to P1 for a "
            f"document-aware re-split before proceeding."
        )

    def scaffold_of(smi):
        m = Chem.MolFromSmiles(smi)
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False) if m else ""

    scaf_by_idx = {i: scaffold_of(src.iloc[i]["canonical_smiles"]) for i in train_kept}
    groups = {}
    for i in train_kept:
        groups.setdefault(scaf_by_idx[i], []).append(i)
    groups_sorted = sorted(groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))

    folds = {str(i): [] for i in range(1, 6)}
    for i, (_, idxs) in enumerate(groups_sorted):
        folds[str((i % 5) + 1)].extend(idxs)

    def n_actives(idxs):
        return int((src.iloc[idxs]["activity_class"] == "active").sum())

    sp["indices"]["train"] = train_kept
    sp["n_compounds"] = {
        "train":      len(train_kept),
        "validation": len(sp["indices"]["validation"]),
        "test":       len(sp["indices"]["test"]),
    }
    sp["n_actives"] = {
        "train":      n_actives(train_kept),
        "validation": n_actives(sp["indices"]["validation"]),
        "test":       n_actives(sp["indices"]["test"]),
    }
    sp["document_leakage_removed_from_train"] = len(train_dropped)
    sp["document_leakage_removed_indices"]    = train_dropped
    cv["folds"] = folds
    print("Document-leakage exclusion applied.")
else:
    print(f"Document-leakage exclusion already applied "
          f"({sp['document_leakage_removed_from_train']} compounds removed). Skipping.")

# ---------------------------------------------------------------
# 2. Rename data/curated_chembl.csv -> eligible-only
# ---------------------------------------------------------------
full_path     = DATA_DIR / "curated_chembl.csv"
eligible_path = DATA_DIR / "curated_chembl_eligible.csv"
archive_path  = DATA_DIR / "curated_chembl_with_ambiguous.csv"

if eligible_path.exists():
    full = pd.read_csv(full_path)
    assert (full["activity_class"] == "ambiguous").any(), \
        "curated_chembl.csv is already eligible-only; skipping rename"
    shutil.copy(full_path, archive_path)
    shutil.copy(eligible_path, full_path)
    eligible_path.unlink()
    sp["curated_csv"] = "data/curated_chembl.csv"
    print(f"Renamed: curated_chembl.csv -> eligible-only ({len(pd.read_csv(full_path))} rows); "
          f"archived full set -> curated_chembl_with_ambiguous.csv")
else:
    print("Eligible file already merged into curated_chembl.csv. Skipping rename.")

# ---------------------------------------------------------------
# 3. Persist split files
# ---------------------------------------------------------------
SP_PATH.write_text(json.dumps(sp, indent=2))
CV_PATH.write_text(json.dumps(cv, indent=2))

# ---------------------------------------------------------------
# 4. Verification
# ---------------------------------------------------------------
print("\n=== VERIFICATION ===")
sp   = json.loads(SP_PATH.read_text())
src2 = pd.read_csv(ROOT / sp["curated_csv"]).reset_index(drop=True)

def docs(idxs):
    return set(src2.iloc[idxs]["document_chembl_id"].dropna())

td, vd, sd = docs(sp["indices"]["train"]), docs(sp["indices"]["validation"]), docs(sp["indices"]["test"])
print(f"curated_csv         : {sp['curated_csv']}")
print(f"train / val / test  : {len(sp['indices']['train'])} / "
      f"{len(sp['indices']['validation'])} / {len(sp['indices']['test'])}")
print(f"train actives       : {sp['n_actives']['train']} / "
      f"{sp['n_actives']['validation']} / {sp['n_actives']['test']} (train/val/test)")
print(f"train ∩ test docs   : {len(td & sd)}  (target 0)")
print(f"train ∩ val  docs   : {len(td & vd)}  (target 0)")
print(f"val   ∩ test docs   : {len(vd & sd)}  (protocol silent — noted)")

for name in ["train", "validation", "test"]:
    cls = src2.iloc[sp["indices"][name]]["activity_class"].value_counts().to_dict()
    assert "ambiguous" not in cls, f"ambiguous found in {name}: {cls}"

print("\nOK")