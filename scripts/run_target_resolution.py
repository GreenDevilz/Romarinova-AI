import os
import sqlite3
from pathlib import Path

import pandas as pd

# Resolve paths relative to the repository root, regardless of CWD
BASE_DIR = Path(__file__).resolve().parent.parent

# Database location — override with CHEMBL_DB_PATH if the dump lives elsewhere
DB_PATH = Path(os.environ.get("CHEMBL_DB_PATH", BASE_DIR / "chembl_37.db")).expanduser()
OUT_CSV = BASE_DIR / "results" / "target_resolution.csv"

if not DB_PATH.exists():
    raise FileNotFoundError(
        f"ChEMBL SQLite database not found at {DB_PATH}.\n"
        "Download the ChEMBL 37 SQLite dump and place it at the repo root as "
        "'chembl_37.db', or set the CHEMBL_DB_PATH environment variable to its location."
    )

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

QUERY = """
SELECT td.tid, td.pref_name, td.organism, td.tax_id,
       COUNT(a.activity_id) AS n_activities
FROM target_dictionary td
LEFT JOIN assays ass     ON ass.tid = td.tid
LEFT JOIN activities a   ON a.assay_id = ass.assay_id
WHERE td.organism = 'Homo sapiens'
  AND td.target_type = 'SINGLE PROTEIN'
  AND (td.pref_name LIKE '%glucosidase%' OR td.pref_name LIKE '%maltase%')
GROUP BY td.tid, td.pref_name, td.organism, td.tax_id
ORDER BY n_activities DESC;
"""

con = sqlite3.connect(DB_PATH)
try:
    df = pd.read_sql_query(QUERY, con)
finally:
    con.close()

df.to_csv(OUT_CSV, index=False)
print(df.to_string(index=False))
print(f"\nWrote {len(df)} rows to {OUT_CSV}")