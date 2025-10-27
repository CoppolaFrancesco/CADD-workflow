#!/usr/bin/env python3
"""
Extract best affinity scores from AutoDock Vina output files and merge with ligand data.
"""

import os
import re
import pandas as pd
from pathlib import Path

def extract_best_affinity(file_path):
    """
    Extract the best (first) affinity value from a vina-score.txt file.
    
    Args:
        file_path: Path to the vina-score.txt file
        
    Returns:
        float: Best affinity value, or None if not found
    """
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Look for the affinity table and extract the first value
        # Pattern matches lines like "   1       -9.723          0          0"
        pattern = r'^\s*1\s+(-?\d+\.\d+)\s+'
        match = re.search(pattern, content, re.MULTILINE)
        
        if match:
            return float(match.group(1))
        else:
            print(f"Warning: Could not find affinity value in {file_path}")
            return None
            
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def process_vina_scores(poses_dir='poses', ligands_csv='ligands/list.csv', output_csv='ligands/list_with_affinities.csv'):
    """
    Process all vina-score.txt files and merge with ligand CSV.
    
    Args:
        poses_dir: Directory containing the vina-score.txt files
        ligands_csv: Path to the input CSV file
        output_csv: Path to save the output CSV file
    """
    
    # Read the ligands CSV
    print(f"Reading ligands data from: {ligands_csv}")
    df = pd.read_csv(ligands_csv)
    
    # Check if 'id-num' column exists (second column should be the ID)
    if 'id-num' not in df.columns:
        print("Warning: 'id-num' column not found. Using second column as ID.")
        id_column = df.columns[1]
    else:
        id_column = 'id-num'
    
    print(f"Using '{id_column}' as the matching column")
    
    # Find all vina-score.txt files
    poses_path = Path(poses_dir)
    
    # Check if directory exists
    if not poses_path.exists():
        print(f"Error: Directory '{poses_dir}' not found!")
        print(f"Current directory: {os.getcwd()}")
        return None
    
    # Try different patterns
    score_files = list(poses_path.glob('*-vina-score.txt'))
    
    # Debug: show what's in the directory
    if len(score_files) == 0:
        all_files = list(poses_path.glob('*'))
        print(f"\nDebug: Found {len(all_files)} files in {poses_dir}/")
        print(f"First 10 files: {[f.name for f in all_files[:10]]}")
        
        # Try finding files in current directory instead
        current_path = Path('.')
        score_files = list(current_path.glob('*-vina-score.txt'))
        if len(score_files) > 0:
            print(f"\nFound {len(score_files)} vina-score.txt files in current directory")
            poses_path = current_path
    
    print(f"Found {len(score_files)} vina-score.txt files")
    
    # Extract affinities
    affinities = {}
    for score_file in score_files:
        # Extract the number from filename (e.g., "0" from "0-vina-score.txt")
        match = re.match(r'(\d+)-vina-score\.txt', score_file.name)
        if match:
            file_num = int(match.group(1))
            affinity = extract_best_affinity(score_file)
            if affinity is not None:
                affinities[file_num] = affinity
    
    print(f"Successfully extracted {len(affinities)} affinity values")
    
    # Create a new column for affinities
    df['vina_affinity'] = df[id_column].map(affinities)
    
    # Check for missing values
    missing_count = df['vina_affinity'].isna().sum()
    if missing_count > 0:
        print(f"Warning: {missing_count} ligands don't have affinity scores")
        missing_ids = df[df['vina_affinity'].isna()][id_column].tolist()
        print(f"Missing IDs: {missing_ids[:10]}..." if len(missing_ids) > 10 else f"Missing IDs: {missing_ids}")
    
    # Save the merged data
    df.to_csv(output_csv, index=False)
    print(f"\nResults saved to: {output_csv}")
    
    # Print summary statistics
    if not df['vina_affinity'].isna().all():
        print("\n" + "="*60)
        print("AFFINITY STATISTICS")
        print("="*60)
        print(f"Total ligands: {len(df)}")
        print(f"Ligands with scores: {len(df) - missing_count}")
        print(f"Best affinity: {df['vina_affinity'].min():.3f} kcal/mol")
        print(f"Worst affinity: {df['vina_affinity'].max():.3f} kcal/mol")
        print(f"Mean affinity: {df['vina_affinity'].mean():.3f} kcal/mol")
        print(f"Median affinity: {df['vina_affinity'].median():.3f} kcal/mol")
        print("="*60)
        
        # Show top 10 best binders
        print("\nTop 10 Best Binders:")
        top_10 = df.nsmallest(10, 'vina_affinity')[[id_column, 'vina_affinity']]
        print(top_10.to_string(index=False))
    
    return df

if __name__ == "__main__":
    # Since you're running from the poses directory, adjust paths
    import sys
    
    # Check if running from poses directory
    current_dir = os.path.basename(os.getcwd())
    
    if current_dir == 'poses':
        print("Running from poses directory")
        df = process_vina_scores(
            poses_dir='.',  # Current directory
            ligands_csv='../ligands/list.csv',
            output_csv='list_with_affinities.csv'  # Output in current (poses) directory
        )
    else:
        print("Running from Autodock-Vina directory")
        df = process_vina_scores(
            poses_dir='poses',
            ligands_csv='ligands/list.csv',
            output_csv='poses/list_with_affinities.csv'
        )
    
    if df is not None:
        print("\n✓ Processing complete!")
    else:
        print("\n✗ Processing failed!")
