import json
import sys
from pathlib import Path

def validate():
    BASE_DIR = Path(__file__).resolve().parent.parent
    contract_path = BASE_DIR / "web" / "data_contract.json"
    
    print(f"Validating {contract_path}...")
    
    if not contract_path.exists():
        print("FAIL: data_contract.json does not exist.")
        sys.exit(1)
        
    try:
        with open(contract_path, "r") as f:
            data = json.load(f)
    except json.JSONDecodeError:
        print("FAIL: data_contract.json is not valid JSON.")
        sys.exit(1)
        
    required_keys = ["metadata", "shortlist_pairs", "compound_level_view", "screened_compounds", "diagnostics"]
    for key in required_keys:
        if key not in data:
            print(f"FAIL: Missing top-level key '{key}'")
            sys.exit(1)
            
    required_meta = ["N_star", "U_star", "known_active_count", "high_ml_ad_extrapolation_count", "ad_extrapolation_count", "ad_missing_count"]
    for key in required_meta:
        if key not in data["metadata"]:
            print(f"FAIL: Missing metadata key '{key}'")
            sys.exit(1)

    if len(data["shortlist_pairs"]) == 0:
        print("FAIL: shortlist_pairs array is empty.")
        sys.exit(1)

    print("SUCCESS: data_contract.json is valid and contains all required fields.")
    sys.exit(0)

if __name__ == "__main__":
    validate()