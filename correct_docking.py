import pandas as pd
import numpy as np
from pathlib import Path

df = pd.read_csv("results/docking_results.csv")

# Apply the new transform
df["docking_score"] = np.clip((-df["docking_affinity_kcal_per_mol"] - 2.0) / 4.0, 0.0, 1.0)

# For failed or positive affinities, force score to 0
df.loc[df["docking_missing"] == True, "docking_score"] = 0.0
df.loc[df["docking_affinity_kcal_per_mol"] > 0, "docking_score"] = 0.0

df.to_csv("results/docking_results.csv", index=False)
print(df[["compound_id", "docking_affinity_kcal_per_mol", "docking_score"]].sort_values("docking_score", ascending=False).head(10).to_string(index=False))