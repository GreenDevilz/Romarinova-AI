"""
Recompute and persist the AD threshold.

The 02_stage1_training notebook computes ad_threshold as a local variable
but does not write it to model_config.json. This script recomputes the
value (95th percentile of training LOO mean Jaccard distance, k=5) from
the locked split and writes it into model_config.json so reviewers can
audit the AD decision.

Run from anywhere:
    python scripts/recover_ad_threshold.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[1]
sp   = json.loads((ROOT / "splits" / "split_indices.json").read_text())
src  = pd.read_csv(ROOT / sp["curated_csv"]).reset_index(drop=True)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
train_idxs = sp["indices"]["train"]
train_fps = [
    gen.GetFingerprint(Chem.MolFromSmiles(src.iloc[i]["canonical_smiles"]))
    for i in train_idxs
]

dists = []
for i, fp in enumerate(train_fps):
    sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
    sims.pop(i)
    sims.sort(reverse=True)
    dists.append(float(np.mean([1 - s for s in sims[:5]])))

threshold = float(np.percentile(dists, 95))

mc_path = ROOT / "models" / "model_config.json"
mc = json.loads(mc_path.read_text()) if mc_path.exists() else {}
mc["ad_threshold"] = threshold
mc["ad_threshold_method"] = (
    "95th percentile of training LOO mean Jaccard distance "
    "(k=5 nearest neighbors, deduplicated training set)"
)
mc["ad_threshold_n_train"] = len(train_idxs)
mc_path.write_text(json.dumps(mc, indent=2))

print(f"AD threshold: {threshold:.4f}")
print(f"  training set size: {len(train_idxs)}")
print(f"  saved to {mc_path}")

# Sanity: report how many rosemary compounds sit near the threshold
ad_path = ROOT / "results" / "ad_results.csv"
if ad_path.exists():
    ad = pd.read_csv(ad_path)
    print(f"  rosemary ad_distance range: "
          f"{ad['ad_distance'].min():.4f} – {ad['ad_distance'].max():.4f}")
    print(f"  compounds within AD: {(ad['ad_flag'] == 'within_AD').sum()} / {len(ad)}")