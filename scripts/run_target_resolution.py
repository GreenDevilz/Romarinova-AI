import sqlite3
from pathlib import Path
import pandas as pd

DB_PATH = Path(r"C:\Users\aymen\chembl_db\chembl_37\chembl_37_sqlite\chembl_37.db")
OUT_CSV = Path("results/target_resolution.csv")
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
df = pd.read_sql_query(QUERY, con)
con.close()

df.to_csv(OUT_CSV, index=False)
print(df.to_string(index=False))
print(f"\nWrote {len(df)} rows to {OUT_CSV}")