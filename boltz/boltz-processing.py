import csv
import os
import subprocess
from pathlib import Path

# Protein sequence
PROTEIN_SEQUENCE = "GPLGSMENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL"

def create_yaml_content(smile):
    """Create YAML content with the given SMILES string"""
    return f"""version: 1  # Optional, defaults to 1
sequences:
  - protein:
      id: A
      sequence: {PROTEIN_SEQUENCE}
  - ligand:
      id: B
      smiles: '{smile}'
properties:
  - affinity:
      binder: B
constraints:
  - pocket:
      binder: B
      contacts: [ [ A, 83 ], [ A, 134 ] ]
"""

def process_ligands(csv_path, output_dir):
    """Process all ligands from the CSV file"""
    
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Read the CSV file
    csv_file = Path(csv_path)
    
    if not csv_file.exists():
        print(f"Error: {csv_file} not found!")
        return
    
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            id_num = row['id-num'].strip()
            smile = row['smiles'].strip()
            
            print(f"\nProcessing ligand {id_num}...")
            
            # Create YAML file
            yaml_path = output_path / f"{id_num}.yaml"
            with open(yaml_path, 'w') as yaml_file:
                yaml_file.write(create_yaml_content(smile))
            
            print(f"Created {yaml_path}")
            
            # Run Boltz prediction
            cmd = [
                "boltz", "predict", str(yaml_path),
                "--use_msa_server",
		"--out_dir", "boltz-results",
                "--recycling_steps", "1",
                "--sampling_steps", "50",
                "--diffusion_samples", "3",
                "--step_scale", "1.2"
            ]
            
            print(f"Running: {' '.join(cmd)}")
            
            try:
                result = subprocess.run(
                    cmd,
                    check=True,
                    capture_output=True,
                    text=True
                )
                print(f"Success! Output:\n{result.stdout}")
            except subprocess.CalledProcessError as e:
                print(f"Error running Boltz for {id_num}:")
                print(f"Return code: {e.returncode}")
                print(f"Stdout: {e.stdout}")
                print(f"Stderr: {e.stderr}")
            except FileNotFoundError:
                print("Error: 'boltz' command not found. Make sure Boltz is installed and in your PATH.")
                return

if __name__ == "__main__":
    # ===================================================================
    # CONFIGURATION - Specify your file paths here
    # ===================================================================
    
    # Path to the input CSV file containing ligands
    CSV_FILE = "../Autodock-Vina/ligands/list.csv"
    
    # Path to the output directory where YAML files will be created
    OUTPUT_DIR = "boltz-configurations-files"
    
    # ===================================================================
    
    process_ligands(CSV_FILE, OUTPUT_DIR)
