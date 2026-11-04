"""
Report scaffold overlap across splits under both Bemis-Murcko definitions.

The split is document-aware (protocol §5.4), which does not guarantee
scaffold-held-out. This quantifies actual cross-split scaffold sharing
so the generalization claim can be stated accurately.

Run from anywhere:
    python scripts/check_scaffold_leakage.py
"""
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]
sp = json.loads((ROOT / "splits" / "split_indices.json").read_text())
src = pd.read_csv(ROOT / sp["curated_csv"]).reset_index(drop=True)


def canonical_scaffold(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)


def generic_scaffold(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    scaf = MurckoScaffold.GetScaffoldForMol(m)
    gen = MurckoScaffold.MakeScaffoldGeneric(scaf)
    return Chem.MolToSmiles(gen, canonical=True)


def report(name, func):
    print(f"\n=== {name} ===")
    scaffolds = [func(s) for s in src["canonical_smiles"]]
    by_split = {}
    for split_name, idxs in sp["indices"].items():
        by_split[split_name] = {scaffolds[i] for i in idxs} - {""}

    total = 0
    for a, b in [("train", "validation"), ("train", "test"), ("validation", "test")]:
        shared = by_split[a] & by_split[b]
        total += len(shared)
        print(f"  {a:12s} ∩ {b:12s}: {len(shared)} shared scaffold(s)")
        for s in sorted(shared)[:5]:
            print(f"    {s}")
        if len(shared) > 5:
            print(f"    ... and {len(shared) - 5} more")
    print(f"  Total: {total}")
    return total


c = report("Canonical Bemis-Murcko (current)", canonical_scaffold)
g = report("Generic Bemis-Murcko", generic_scaffold)

print(f"\n--- Summary ---")
print(f"Canonical overlap: {c}")
print(f"Generic overlap:   {g}")
print(f"Declared in split_indices.json: {sp.get('scaffold_definition')}")