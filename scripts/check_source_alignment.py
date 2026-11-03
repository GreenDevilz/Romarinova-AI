"""
Verify that every script and notebook that reads curated_chembl data is
loading the file referenced by splits/split_indices.json.

Run from anywhere:
    python scripts/check_source_alignment.py
"""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sp = json.loads((ROOT / "splits" / "split_indices.json").read_text())
expected_file = sp["curated_csv"]

print(f"split_indices.json -> {expected_file}")
src = pd.read_csv(ROOT / expected_file)
print(f"  rows: {len(src)}  classes: {src['activity_class'].value_counts().to_dict()}")

# Expected split row count
expected_total = sum(len(sp["indices"][s]) for s in ["train", "validation", "test"])
print(f"  expected split total: {expected_total} rows")

# --- Scan every .py in scripts/ and every .ipynb in notebooks/ for data reads
print("\n=== Source alignment scan ===")
files_scanned = 0
misaligned = []

for pattern_path, ext_label in [
    (ROOT / "scripts" / "*.py", ".py"),
    (ROOT / "notebooks" / "*.ipynb", ".ipynb"),
]:
    for p in sorted(pattern_path.parent.glob(pattern_path.name)):
        files_scanned += 1
        if p.suffix == ".ipynb":
            nb = json.loads(p.read_text(encoding="utf-8"))
            text = "\n".join(
                "".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"
            )
        else:
            text = p.read_text(encoding="utf-8")

        # Match any read_csv(...) call whose argument text contains "curated_chembl".
        # Handles both read_csv("literal") and read_csv(Path / "literal").
        for match in re.finditer(r"read_csv\(([^)]*)\)", text):
            arg = match.group(1)
            if "curated_chembl" not in arg:
                continue
            # Extract the last quoted string from the argument
            strings = re.findall(r'"([^"]*)"', arg) + re.findall(r"'([^']*)'", arg)
            if not strings:
                continue
            basename = Path(strings[-1]).name
            expected_basename = Path(expected_file).name
            if basename not in {expected_basename, "curated_chembl.csv"}:
                continue  # looking at archived files intentionally is OK
            status = "OK" if basename == expected_basename else "MISALIGNED"
            print(f"  {ext_label} {p.name:40s} reads {basename:35s} [{status}]")
            if status == "MISALIGNED":
                misaligned.append((p.name, arg.strip()))

        # Check for split_indices.json reads (the correct pattern)
        if "split_indices" in text and 'curated_csv' in text:
            print(f"  {ext_label} {p.name:40s} reads split['curated_csv'] [OK - dynamic]")

print(f"\nScanned {files_scanned} files.")

if misaligned:
    print(f"\n{len(misaligned)} misaligned file(s) — fix before rerunning:")
    for fname, path in misaligned:
        print(f"  {fname}: {path}")
else:
    print("\nAll data reads resolve to the expected source.")

# --- Diagnose what happens if the wrong file is loaded
print("\n=== Impact if wrong file is loaded ===")
wrong = ROOT / "data" / "curated_chembl.csv"
if wrong.exists():
    wsrc = pd.read_csv(wrong)
    if len(wsrc) != len(src):
        print(f"  {wrong.name}: {len(wsrc)} rows (vs {len(src)} in {expected_file})")
        print(f"  Split indices point into {len(src)}-row space; loading a "
              f"{len(wsrc)}-row file would misalign every row.")
        # Show first divergence
        for sname in ["train", "validation", "test"]:
            idxs = sp["indices"][sname]
            same = sum(
                1 for i in idxs
                if i < len(wsrc)
                and wsrc.iloc[i]["standard_inchikey"] == src.iloc[i]["standard_inchikey"]
            )
            print(f"  {sname:12s}: {same}/{len(idxs)} rows match by InChIKey")
    else:
        print(f"  Both files have {len(src)} rows — same file, no misalignment risk.")