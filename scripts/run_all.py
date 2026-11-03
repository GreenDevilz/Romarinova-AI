import subprocess
import sys

def run_step(script_path, description):
    print(f"\n{'='*40}")
    print(f"Executing: {description}")
    print(f"Command: python {script_path}")
    print(f"{'='*40}")
    
    result = subprocess.run([sys.executable, script_path], text=True)
    
    if result.returncode != 0:
        print(f"\nFAILED: {script_path} exited with error code {result.returncode}.")
        sys.exit(1)
    print(f"SUCCESS: {script_path} completed.")

if __name__ == "__main__":
    print("Starting Romarinova AI End-to-End Pipeline...")
    
    # Note: chembl_query.py and target_resolution.py require the massive DB, 
    # so we start the automated run from the artifacts step.
    
    run_step("scripts/evaluation.py", "Model Evaluation Table")
    run_step("scripts/sensitivity.py", "Sensitivity Analysis (Ligand-Only vs Locked)")
    run_step("scripts/build_artifacts.py", "Generate Frontend Payload & Figures")
    run_step("scripts/validate_data_contract.py", "Validate Data Contract for Frontend")
    
    print("\nPipeline complete. Ready for frontend deployment.")