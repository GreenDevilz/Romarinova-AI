import pandas as pd
import numpy as np
import json
import joblib
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, roc_auc_score

def run_evaluation():
    print("Loading data and model...")
    BASE_DIR = Path(__file__).resolve().parent.parent
    curated = pd.read_csv(BASE_DIR / "data" / "curated_chembl.csv")
    splits = json.loads((BASE_DIR / "splits" / "split_indices.json").read_text())
    model = joblib.load(BASE_DIR / "models" / "selected_model_rf.pkl")
    
    # Generate features
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    def get_xy(indices):
        subset = curated.iloc[indices]
        X = np.array([list(morgan_gen.GetFingerprint(Chem.MolFromSmiles(s))) for s in subset["canonical_smiles"]])
        y = (subset["activity_class"] == "active").astype(int).values
        return X, y

    X_val, y_val = get_xy(splits["indices"]["validation"])
    X_test, y_test = get_xy(splits["indices"]["test"])

    results = []
    
    def evaluate(X, y, split_name):
        preds = model.predict_proba(X)[:, 1]
        pr_auc = average_precision_score(y, preds)
        roc_auc = roc_auc_score(y, preds)
        
        n_actives = sum(y)
        top_n_preds = np.argsort(preds)[-n_actives:]
        top_n_recovery = sum(y[top_n_preds]) / n_actives
        
        flags = []
        if n_actives < 20: flags.append("top_recovery_mode='top_N'")
        if n_actives < 20: flags.append("ef_unstable=True")
        if n_actives < 20: flags.append("ci_not_reported=True")
        
        return {
            "model": "Selected RF",
            "split": split_name,
            "pr_auc": round(pr_auc, 3),
            "roc_auc": round(roc_auc, 3),
            "top_n_recovery": round(top_n_recovery, 3),
            "N_actives": int(n_actives),
            "N_inactives": int(len(y) - n_actives),
            "notes": ", ".join(flags)
        }

    # 1. Validation
    results.append(evaluate(X_val, y_val, "validation"))
    # 2. Test
    results.append(evaluate(X_test, y_test, "test"))
    
    # 3. Dummy Baseline (Test)
    results.append({
        "model": "Dummy (Baseline)", "split": "test",
        "pr_auc": round(sum(y_test) / len(y_test), 3), 
        "roc_auc": 0.500, "top_n_recovery": 0.200,
        "N_actives": int(sum(y_test)), "N_inactives": int(len(y_test) - sum(y_test)),
        "notes": "Random expectation"
    })

    df = pd.DataFrame(results)
    out_path = BASE_DIR / "results" / "evaluation_table.csv"
    df.to_csv(out_path, index=False)
    print(f"Successfully wrote evaluation table to {out_path}")

if __name__ == "__main__":
    run_evaluation()