import os
import sqlite3
from pathlib import Path

import pandas as pd

# Resolve paths relative to the repository root, regardless of CWD
BASE_DIR = Path(__file__).resolve().parent.parent

# Database location — override with CHEMBL_DB_PATH if the dump lives elsewhere
DB_PATH = Path(os.environ.get("CHEMBL_DB_PATH", BASE_DIR / "chembl_37.db")).expanduser()
SQL_PATH = BASE_DIR / "scripts" / "chembl_query.sql"
OUT_PATH = BASE_DIR / "data" / "raw_chembl_tid1.csv"

if not DB_PATH.exists():
    raise FileNotFoundError(
        f"ChEMBL SQLite database not found at {DB_PATH}.\n"
        "Download the ChEMBL 37 SQLite dump and place it at the repo root as "
        "'chembl_37.db', or set the CHEMBL_DB_PATH environment variable to its location."
    )
if not SQL_PATH.exists():
    raise FileNotFoundError(f"SQL query file not found at {SQL_PATH}.")

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

sql = SQL_PATH.read_text(encoding="utf-8")
con = sqlite3.connect(DB_PATH)
try:
    df = pd.read_sql_query(sql, con)
finally:
    con.close()

df.to_csv(OUT_PATH, index=False)
print(f"Wrote {len(df)} rows × {len(df.columns)} cols to {OUT_PATH}")
print(df["activity_type"].value_counts(dropna=False).to_string())
print("\nMissing InChIKey rows:", df["standard_inchikey"].isna().sum())
print("Missing SMILES rows:", df["canonical_smiles"].isna().sum())