import sqlite3
from pathlib import Path
import pandas as pd

DB_PATH  = Path(r"C:\Users\aymen\chembl_db\chembl_37\chembl_37_sqlite\chembl_37.db")
SQL_PATH = Path("scripts/chembl_query.sql")
OUT_PATH = Path("data/raw_chembl_tid1.csv")
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

sql = SQL_PATH.read_text(encoding="utf-8")
con = sqlite3.connect(DB_PATH)
df = pd.read_sql_query(sql, con)
con.close()

df.to_csv(OUT_PATH, index=False)
print(f"Wrote {len(df)} rows × {len(df.columns)} cols to {OUT_PATH}")
print(df["activity_type"].value_counts(dropna=False).to_string())
print("\nMissing InChIKey rows:", df["standard_inchikey"].isna().sum())
print("Missing SMILES rows:", df["canonical_smiles"].isna().sum())