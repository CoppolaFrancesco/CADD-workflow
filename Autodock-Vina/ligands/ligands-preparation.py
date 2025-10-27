import csv
import subprocess
import os
import glob

def add_id_column_to_cheese_file(target_folder):
    """
    Find CSV files starting with 'cheese', add an 'id-num' column starting from 0,
    and save the result as 'list.csv'.
    
    Args:
        target_folder (str): Path to the folder to search for cheese*.csv files
        
    Returns:
        str: Name of the output file ('list.csv'), or None if no file found
    """
    # Find all CSV files starting with 'cheese'
    search_pattern = os.path.join(target_folder, 'cheese*.csv')
    cheese_files = glob.glob(search_pattern)
    
    if not cheese_files:
        print("No CSV files starting with 'cheese' found in the folder.")
        return None
    
    if len(cheese_files) > 1:
        print(f"Warning: Found {len(cheese_files)} files matching 'cheese*.csv'. Using the first one: {cheese_files[0]}")
    
    cheese_file = cheese_files[0]
    output_file = os.path.join(target_folder, 'list.csv')
    
    print(f"--- Processing cheese file: {cheese_file} ---")
    print(f"--- Output will be saved to: {output_file} ---")
    
    # Read the original file
    rows = []
    try:
        with open(cheese_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            rows = list(reader)
    except Exception as e:
        print(f"Error reading {cheese_file}: {e}")
        return None
    
    if not rows:
        print(f"File {cheese_file} is empty.")
        return None
    
    # Add 'id-num' as the second column (after SMILES, before any other columns)
    new_rows = []
    
    # Process header
    header = rows[0] if rows else []
    if header:
        new_header = [header[0], 'id-num'] + header[1:]
        new_rows.append(new_header)
    
    # Process data rows - add id starting from 0
    for idx, row in enumerate(rows[1:], start=0):
        if row:  # Skip completely empty rows
            new_row = [row[0], str(idx)] + row[1:]
            new_rows.append(new_row)
    
    # Write to list.csv
    try:
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(new_rows)
        print(f"✓ Created {output_file} with 'id-num' column added")
        print(f"  Total data rows: {len(new_rows) - 1}\n")
        return 'list.csv'
    except Exception as e:
        print(f"Error writing to {output_file}: {e}")
        return None


def process_smiles_file(target_folder, target_file):
    """
    Process a CSV file containing SMILES strings through a two-step ligand preparation pipeline.
    
    Step 1: scrub.py - Converts SMILES to SDF format
    Step 2: mk_prepare_ligand.py - Converts SDF to PDBQT format
    
    Args:
        target_folder (str): Path to the folder containing the CSV file
        target_file (str): Name of the CSV file to process
    """
    
    csv_file_path = os.path.join(target_folder, target_file)
    
    if not os.path.exists(csv_file_path):
        print(f"ERROR: CSV file not found at '{csv_file_path}'")
        return
    
    print(f"--- Starting Ligand Preparation Pipeline for: {csv_file_path} ---\n")
    
    # Track statistics
    processed_count = 0
    successful_conversions = 0
    failed_conversions = []
    consecutive_empty_rows = 0
    MAX_CONSECUTIVE_EMPTY = 5  # Stop if we hit this many empty rows in a row
    all_processed_ids = []  # Track IDs we actually attempted to process
    
    try:
        with open(csv_file_path, 'r', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            
            # Read and skip the header row
            header = next(reader, None)
            if header:
                print(f"Header row: {header}\n")
            
            # Process each subsequent row until the file ends or we hit too many empty rows
            for row_num, row in enumerate(reader, start=2):
                # Skip empty or malformed rows
                if not row or len(row) < 2:
                    consecutive_empty_rows += 1
                    if consecutive_empty_rows >= MAX_CONSECUTIVE_EMPTY:
                        print(f"\nEncountered {MAX_CONSECUTIVE_EMPTY} consecutive empty rows. Stopping processing.")
                        break
                    continue
                
                smiles_string = row[0].strip()
                ligand_id = row[1].strip()
                
                # Skip rows where either SMILES or ligand_id is empty
                if not smiles_string or not ligand_id:
                    consecutive_empty_rows += 1
                    if consecutive_empty_rows >= MAX_CONSECUTIVE_EMPTY:
                        print(f"\nEncountered {MAX_CONSECUTIVE_EMPTY} consecutive empty rows. Stopping processing.")
                        break
                    continue
                
                # Reset counter when we find valid data
                consecutive_empty_rows = 0
                processed_count += 1
                all_processed_ids.append(ligand_id)  # Track this ID
                
                # Progress indicator - just show count
                print(f"[Ligand {processed_count}] Processing ID: {ligand_id}")
                print(f"  SMILES: {smiles_string}")
                
                # Define output filenames
                sdf_output_filename = f"{ligand_id}-prepared.sdf"
                pdbqt_output_filename = f"{ligand_id}-prepared.pdbqt"
                
                scrub_success = False
                
                # --- STEP 1: Run scrub.py ---
                try:
                    print(f"  > Step 1: Running scrub.py...")
                    command_scrub = ['scrub.py', smiles_string, '-o', sdf_output_filename]
                    
                    result_scrub = subprocess.run(
                        command_scrub,
                        capture_output=True, text=True, check=False
                    )

                    if result_scrub.returncode == 0:
                        print(f"  > scrub.py SUCCESS. SDF file created: '{sdf_output_filename}'")
                        scrub_success = True
                    else: 
                        print(f"  > scrub.py FAILED (Exit Code: {result_scrub.returncode}). Skipping Step 2.")
                        print(f"  > Stderr from scrub.py: {result_scrub.stderr.strip()}")
                        failed_conversions.append((ligand_id, "scrub.py failed"))
                        
                except FileNotFoundError:
                    print(f"  > CRITICAL ERROR: 'scrub.py' not found. Check your PATH.")
                    failed_conversions.append((ligand_id, "scrub.py not found"))
                    break
                except Exception as e:
                    print(f"  > An unexpected error occurred during scrub.py execution: {e}")
                    failed_conversions.append((ligand_id, f"scrub.py error: {e}"))
       
                # --- STEP 2: Run mk_prepare_ligand.py if scrub.py was successful ---
                if scrub_success:
                    try:
                        print(f"  > Step 2: Running mk_prepare_ligand.py...")
                        command_mk = ['mk_prepare_ligand.py', '-i', sdf_output_filename, '-o', pdbqt_output_filename]
                        
                        result_mk = subprocess.run(
                            command_mk,
                            capture_output=True, text=True, check=False
                        )
        
                        if result_mk.returncode == 0:
                            print(f"  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '{pdbqt_output_filename}'")
                            successful_conversions += 1
                        else:
                            print(f"  > mk_prepare_ligand.py FAILED (Exit Code: {result_mk.returncode}).")
                            print(f"  > Stderr from mk_prepare_ligand.py: {result_mk.stderr.strip()}")
                            failed_conversions.append((ligand_id, "mk_prepare_ligand.py failed"))
                            
                    except FileNotFoundError:
                        print(f"  > CRITICAL ERROR: 'mk_prepare_ligand.py' not found. Check your PATH.")
                        failed_conversions.append((ligand_id, "mk_prepare_ligand.py not found"))
                        break
                    except Exception as e:
                        print(f"  > An unexpected error occurred during mk_prepare_ligand.py execution: {e}")
                        failed_conversions.append((ligand_id, f"mk_prepare_ligand.py error: {e}"))
                
                print()  # Blank line between ligands
            
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")

    print("\n" + "="*70)
    print("--- PIPELINE VERIFICATION REPORT ---")
    print("="*70)
    
    # Verify files exist on disk for the ligands we actually processed
    print("\n1. FILE EXISTENCE CHECK:")
    print("-" * 70)
    
    missing_sdf = []
    missing_pdbqt = []
    
    # Check only the IDs we actually processed
    for ligand_id in all_processed_ids:
        sdf_file = f"{ligand_id}-prepared.sdf"
        pdbqt_file = f"{ligand_id}-prepared.pdbqt"
        
        if not os.path.exists(sdf_file):
            missing_sdf.append(ligand_id)
        
        if not os.path.exists(pdbqt_file):
            missing_pdbqt.append(ligand_id)
    
    # Summary statistics
    print(f"\n2. PROCESSING SUMMARY:")
    print("-" * 70)
    print(f"Total ligands processed:          {len(all_processed_ids)}")
    print(f"Successfully converted to PDBQT:  {successful_conversions}")
    print(f"Failed conversions:               {len(failed_conversions)}")
    print(f"Missing SDF files:                {len(missing_sdf)}")
    print(f"Missing PDBQT files:              {len(missing_pdbqt)}")
    
    # Success rate
    if len(all_processed_ids) > 0:
        success_rate = (successful_conversions / len(all_processed_ids)) * 100
        print(f"\nSuccess rate:                     {success_rate:.1f}%")
    
    # Detailed failure report
    if failed_conversions:
        print(f"\n3. FAILED CONVERSIONS DETAILS:")
        print("-" * 70)
        for ligand_id, reason in failed_conversions:
            print(f"  • ID {ligand_id}: {reason}")
    
    # Missing files report
    if missing_sdf:
        print(f"\n4. MISSING SDF FILES:")
        print("-" * 70)
        for ligand_id in missing_sdf:
            print(f"  • {ligand_id}-prepared.sdf")
    
    if missing_pdbqt:
        print(f"\n5. MISSING PDBQT FILES:")
        print("-" * 70)
        for ligand_id in missing_pdbqt:
            print(f"  • {ligand_id}-prepared.pdbqt")
    
    # Final verdict
    print("\n" + "="*70)
    if successful_conversions == len(all_processed_ids) and not missing_pdbqt:
        print("✓ SUCCESS: All ligands were successfully converted!")
    else:
        print("✗ WARNING: Some ligands failed to convert. Review details above.")
    print("="*70)
    
    print(f"\n--- Pipeline finished: Processed {processed_count} ligands ---")
                    
# --- Configuration Section ---
TARGET_FOLDER = '.'
                
if __name__ == '__main__':
    # Step 1: Find and process cheese*.csv file (add id-num column and save as list.csv)
    print("="*70)
    print("STEP 1: Looking for cheese*.csv file to create list.csv")
    print("="*70 + "\n")
    
    result_file = add_id_column_to_cheese_file(TARGET_FOLDER)
    
    if not result_file:
        print("ERROR: No cheese*.csv file found. Cannot proceed with pipeline.")
        print("Please ensure there is a CSV file starting with 'cheese' in the folder.\n")
        exit(1)
    
    # Step 2: Process the list.csv file through the ligand preparation pipeline
    print("="*70)
    print("STEP 2: Ligand Preparation Pipeline")
    print("="*70 + "\n")
    
    process_smiles_file(TARGET_FOLDER, 'list.csv')
