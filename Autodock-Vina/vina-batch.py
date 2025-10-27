#!/usr/bin/env python3
"""
AutoDock Vina Batch Docking Script
Performs docking for multiple ligands and saves results organized in poses directory
"""

import os
import subprocess
from pathlib import Path

def run_vina_docking(ligands_dir='ligands', 
                     receptor_dir='receptor',
                     receptor_name='1H1Q-prepared.pdbqt',
                     config_file='1H1Q-prepared.box.txt',
                     exhaustiveness=100,
                     num_modes=20,
                     poses_dir='poses'):
    """
    Run AutoDock Vina docking for all ligands in the specified directory.
    
    Args:
        ligands_dir: Directory containing ligand PDBQT files
        receptor_dir: Directory containing receptor files
        receptor_name: Name of the receptor PDBQT file
        config_file: Name of the box configuration file
        exhaustiveness: Exhaustiveness parameter for Vina
        num_modes: Number of binding modes to generate
        poses_dir: Directory to save output files and poses
    """
    
    # Setup paths
    receptor_path = os.path.join(receptor_dir, receptor_name)
    config_path = os.path.join(receptor_dir, config_file)
    
    # Verify required files exist
    if not os.path.exists(receptor_path):
        print(f"Error: Receptor file not found: {receptor_path}")
        return
    
    if not os.path.exists(config_path):
        print(f"Error: Config file not found: {config_path}")
        return
    
    if not os.path.exists(ligands_dir):
        print(f"Error: Ligands directory not found: {ligands_dir}")
        return
    
    # Create poses directory if it doesn't exist
    os.makedirs(poses_dir, exist_ok=True)
    
    # Get all PDBQT ligand files and extract ligand numbers
    ligand_files = sorted([f for f in os.listdir(ligands_dir) if f.endswith('-prepared.pdbqt')])
    
    if not ligand_files:
        print(f"Error: No ligand files found in {ligands_dir}")
        return
    
    # Extract ligand info (number, full path)
    ligand_info = []
    for lf in ligand_files:
        try:
            ligand_num = lf.split('-')[0]
            full_path = os.path.join(ligands_dir, lf)
            ligand_info.append((ligand_num, full_path))
        except:
            print(f"Warning: Could not parse ligand number from {lf}")
            continue
    
    print(f"Found {len(ligand_info)} ligand files")
    print("Starting docking process...\n")
    
    # Track results
    successful = []
    failed = []
    
    # Run Vina for each ligand
    for ligand_num, ligand_path in ligand_info:
        # Save output files in poses directory
        output_file = os.path.join(poses_dir, f"{ligand_num}-vina-score.txt")
        output_pose = os.path.join(poses_dir, f"{ligand_num}-vina-out.pdbqt")
                     
        print(f"Processing ligand {ligand_num}...")
       
        # Build the vina command
        cmd = [
            'vina',
            '--receptor', receptor_path,
            '--ligand', ligand_path,
            '--config', config_path,
            f'--exhaustiveness={exhaustiveness}',
            '--out', output_pose,
            '--num_modes', str(num_modes)
        ]
       
        try:
            # Run the command and redirect output to file
            with open(output_file, 'w') as out_f:
                result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE,
                                       text=True, check=True)

            print(f"  ✓ Completed: {output_file}")
            successful.append(ligand_num)
            
        except subprocess.CalledProcessError as e:
            print(f"  ✗ Failed: Ligand {ligand_num}")
            print(f"    Error: {e.stderr}")
            failed.append(ligand_num)
        except FileNotFoundError:
            print(f"  ✗ Error: 'vina' command not found. Make sure AutoDock Vina is installed and in your PATH")
            return

    # Final verification and summary
    print(f"\n{'='*70}")
    print("DOCKING SUMMARY")
    print(f"{'='*70}\n")
    
    # Count expected vs actual output files
    expected_score_files = len(ligand_info)
    expected_pose_files = len(ligand_info)
    
    actual_score_files = len([f for f in os.listdir(poses_dir) if f.endswith('-vina-score.txt')])
    actual_pose_files = len([f for f in os.listdir(poses_dir) if f.endswith('-vina-out.pdbqt')])
    
    print(f"Total ligands processed:     {len(ligand_info)}")
    print(f"Successful dockings:         {len(successful)}")
    print(f"Failed dockings:             {len(failed)}")
    print(f"\nOutput files verification:")
    print(f"  Score files (.txt):        {actual_score_files}/{expected_score_files}")
    print(f"  Pose files (.pdbqt):       {actual_pose_files}/{expected_pose_files}")
    
    success_rate = (len(successful) / len(ligand_info) * 100) if ligand_info else 0
    print(f"\nSuccess rate:                {success_rate:.1f}%")
    
    if failed:
        print(f"\nFailed ligands: {', '.join(failed)}")
    
    print(f"\nAll outputs saved in: {os.path.abspath(poses_dir)}/")
    print(f"{'='*70}")
    
    # Final status
    if len(successful) == len(ligand_info):
        print("\n✓ SUCCESS: All ligands docked successfully!")
    elif len(successful) > 0:
        print(f"\n⚠ WARNING: {len(failed)} ligand(s) failed to dock")
    else:
        print("\n✗ ERROR: All dockings failed")
    
    print(f"{'='*70}\n")

if __name__ == "__main__":
    # Run the docking with default parameters
    run_vina_docking(
        ligands_dir='ligands',
        receptor_dir='receptor',
        receptor_name='1H1Q-prepared.pdbqt',
        config_file='1H1Q-prepared.box.txt',
        exhaustiveness=100,
        num_modes=20,
        poses_dir='poses'
    )
