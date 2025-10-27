import os
import json
import pandas as pd
import numpy as np

def analyze_boltz_results(results_dir, csv_path, output_path):
    """
    Analyze Boltz prediction results and merge with existing CSV
    
    Args:
        results_dir: Directory containing boltz_results_* folders
        csv_path: Path to the original CSV with affinities
        output_path: Path where the new CSV will be saved
    """
    
    # ------------------------------------------------------------ 
    # 1. Load original CSV
    # ------------------------------------------------------------ 
    if not os.path.exists(csv_path):
        print(f"❌ Error: CSV file not found at {csv_path}")
        return
    
    df = pd.read_csv(csv_path)
    print(f"📄 Loaded CSV with {len(df)} rows")
    
    # Dictionary to store boltz data
    boltz_affinities = {}
    avg_pred_values = {}
    avg_prob_binary = {}
    
    # ------------------------------------------------------------ 
    # 2. Loop through all boltz_results_* folders
    # ------------------------------------------------------------ 
    if not os.path.exists(results_dir):
        print(f"❌ Error: Results directory not found at {results_dir}")
        return
    
    processed_count = 0
    missing_count = 0
    
    for folder_name in os.listdir(results_dir):
        if folder_name.startswith("boltz_results_"):
            try:
                # Extract numeric suffix or full ID
                idx = folder_name.replace("boltz_results_", "")
                # Try to convert to int if possible
                try:
                    idx = int(idx)
                except ValueError:
                    # Keep as string if not a simple integer
                    pass
            except ValueError:
                continue
            
            # Build JSON path: boltz-results/boltz_results_X/predictions/X/affinity_X.json
            json_path = os.path.join(results_dir, folder_name, "predictions", str(idx), f"affinity_{idx}.json")
            
            if not os.path.exists(json_path):
                print(f"⚠️  Warning: Missing file {json_path}")
                missing_count += 1
                continue
            
            # ------------------------------------------------------------ 
            # 3. Read JSON and compute averages
            # ------------------------------------------------------------ 
            try:
                with open(json_path, "r") as f:
                    data = json.load(f)
                
                pred_values = [
                    data["affinity_pred_value"],
                    data["affinity_pred_value1"],
                    data["affinity_pred_value2"]
                ]
                
                prob_binary_values = [
                    data["affinity_probability_binary"],
                    data["affinity_probability_binary1"],
                    data["affinity_probability_binary2"]
                ]
                
                mean_pred_value = np.mean(pred_values)
                mean_prob_binary = np.mean(prob_binary_values)
                
                # ------------------------------------------------------------ 
                # 4. Convert mean affinity to kcal/mol
                # ------------------------------------------------------------ 
                boltz_kcalmol = (6 - mean_pred_value) * 1.364
                boltz_affinities[idx] = boltz_kcalmol
                avg_pred_values[idx] = mean_pred_value
                avg_prob_binary[idx] = mean_prob_binary
                
                processed_count += 1
                print(f"✅ Processed {folder_name}: {boltz_kcalmol:.2f} kcal/mol (pred={mean_pred_value:.3f}, prob={mean_prob_binary:.3f})")
                
            except KeyError as e:
                print(f"⚠️  Missing values in {json_path}: {e}")
                missing_count += 1
                continue
            except json.JSONDecodeError as e:
                print(f"⚠️  Error parsing JSON in {json_path}: {e}")
                missing_count += 1
                continue
    
    print(f"\n📊 Summary: Processed {processed_count} results, {missing_count} missing/errors")
    
    # ------------------------------------------------------------ 
    # 5. Merge with CSV (match by ligand ID = 2nd column)
    # ------------------------------------------------------------ 
    id_col = df.columns[1]
    print(f"🔗 Matching on column: {id_col}")
    
    df["boltz_affinity_kcalmol"] = df[id_col].map(boltz_affinities)
    df["avg_affinity_pred_value"] = df[id_col].map(avg_pred_values)
    df["avg_affinity_probability_binary"] = df[id_col].map(avg_prob_binary)
    
    # Count how many matches were found
    matched = df["boltz_affinity_kcalmol"].notna().sum()
    print(f"🎯 Matched {matched}/{len(df)} rows with Boltz results")
    
    # ------------------------------------------------------------ 
    # 6. Save new CSV
    # ------------------------------------------------------------ 
    df.to_csv(output_path, index=False)
    print(f"\n✅ Done! New CSV saved to: {output_path}")
    print(f"📈 Added columns: boltz_affinity_kcalmol, avg_affinity_pred_value, avg_affinity_probability_binary")
    
    # Show some statistics
    if matched > 0:
        print(f"\n📉 Boltz Affinity Statistics:")
        print(f"   Mean: {df['boltz_affinity_kcalmol'].mean():.2f} kcal/mol")
        print(f"   Std:  {df['boltz_affinity_kcalmol'].std():.2f} kcal/mol")
        print(f"   Min:  {df['boltz_affinity_kcalmol'].min():.2f} kcal/mol")
        print(f"   Max:  {df['boltz_affinity_kcalmol'].max():.2f} kcal/mol")

if __name__ == "__main__":
    # ===================================================================
    # CONFIGURATION - Specify your file paths here
    # ===================================================================
    
    # Directory containing boltz_results_* folders
    RESULTS_DIR = "boltz-results"
    
    # Path to the original CSV file with affinities from AutoDock Vina
    INPUT_CSV = "../Autodock-Vina/poses/list_with_affinities.csv"
    
    # Path where the new CSV with Boltz results will be saved
    OUTPUT_CSV = "list_with_affinities_boltz.csv"
    
    # ===================================================================
    
    analyze_boltz_results(RESULTS_DIR, INPUT_CSV, OUTPUT_CSV)
