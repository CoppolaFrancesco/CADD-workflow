# Documentation

## Table of Contents
* [Detailed Documentation](#detailed-documentation)
* [Files Preparation](#files-preparation)
  * [Receptor Preparation](#receptor-preparation)
  * [Ligand Preparation](#ligands-preparation)
* [AutoDock Vina](#autodock-vina)
  * [Automatization Vina](#automatization-vina)
* [Boltz](#boltz)
  * [Automatization Boltz](#automatization-boltz)
* [Additional Properties](#additional-properties)
  * [Sorting](#sorting)
  * [RDKit and DeepChem](#rdkit-and-deepchem)

# Detailed Documentation
In this document, we will execute each step manually, without relying on the automated bash scripts. All necessary files are included in the folders described below. This setup workflow enables me to manage each step individually and easily identify potential problems. However, these steps could be combined into a single file, such as consolidating all the Python scripts. The original idea was to offer several alternatives. For example, instead of using AutoDock Vina, one could run Uni-Dock. Creating a separate block for this alternative during the design of the workflow provides several benefits, enhancing the overall flexibility of the workflow. Another interesting feature of this architecture is that it allows us to easily add a parallelization process in the future, where we can create and process a batch of several ligands on different devices without big modifications.

Every step is logged, and I will report in this documentation the outputs. 

# Files Preparation

## Receptor Preparation

In this demonstration, I will be using CDK2 as the target protein, which is represented by the PDB code 1H1Q. The PDB files typically contain many water molecules (some of which are important!!), cofactors, other proteins, and ligands. For now, for simplicity, we are going to take only the receptor. I suggest two approaches: one is to use [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) or the PDB Reader & Manipulator in [CHARMM-GUI](https://www.charmm-gui.org/?doc=input/pdbreader). I'm currently using the second option, and in Step 1 of CHARMM-GUI, we will be prompted to select the chain we want to model. I selected to model only chain A (segid PROA). Step 2 is particularly important, as it allows us to manipulate the PDB file. We need to carefully consider the protonation states and any missing residues. Once we have completed these steps, we can download the PDB file. I have also included this file in the folder `Autodock-Vina/receptor/1H1Q-receptor-only.pdb`, along with the other generated files.

We need to generate the PDBQT file for our receptor, which is a modified .pdb file that includes partial charges and atom types. Using Meeko (which was already installed in the [introduction](../README.md)), we can generate the PDBQT file of our receptor along with the docking pocket box. I approximated the dimensions of the box based on the crystallographic structure of the complex with the ligand NU6094. This box will define the area for our docking search. 

```
# make sure you are in the vina environment
conda activate vina
mk_prepare_receptor.py -i 1H1Q-receptor-only.pdb -o 1H1Q-prepared -p -v --box_size 20 20 20 --box_center 6.145 44.175 50.827 
```

In my case, I plan to perform a semi-flexible docking where the receptor will remain completely frozen. If we wanted to conduct flexible docking, we would simply add the flag `-f A:83` to the previous command to make residue 83 of chain A flexible, for example. However, since we are using the crystallographic structure of the protein that is already bound to the inhibitor NU6094, I will keep both the pocket and the receptor frozen. Another reason to keep the pocket frozen is that my list of ligands is similar to NU6094. If we were using a wider variety of molecules with less similarity, it would likely be more beneficial to allow the pocket to be flexible. Finally, conducting a thorough literature review is essential to deciding how to proceed in this case. The conformational changes of the protein are also a crucial feature of our docking, especially when dealing with allosteric inhibitors.

The previous command should now generate the following files `1H1Q-prepared.pdbqt`, `1H1Q-prepared.box.txt` and `1H1Q-prepared.box.pdb` and prints on our console:

```
@> 2388 atoms and 1 coordinate set(s) were parsed in 0.01s.

Files written:
  1H1Q-prepared.pdbqt <-- static (i.e., rigid) receptor input file
1H1Q-prepared.box.txt <-- Vina-style box dimension file
1H1Q-prepared.box.pdb <-- PDB file to visualize the grid box
```

## Ligands Preparation

For the dataset preparation, I used the [cheese search](https://cheese.deepmedchem.com/). Specifically, by searching exclusively in the Mcule Full database, using a method known as [Shape: 3D volume overlap](https://chemrxiv.org/engage/chemrxiv/article-details/67250915f9980725cfcd1f6f). In the dataset, I aimed to include a few analogs of the ligand (NU6094) as well as a range of diverse compounds with different structural similarities and numbers of rotatable bonds. I selected approximately 300 molecules and exported them as a CSV file. The CSV file is structured so that the first column contains the SMILES representation. We need to review each smile, add the hydrogens using Molscrub, save the result as a .sdf file, and then generate .pdbqt files with Meeko. 

These are the two ideally steps that one should run:

```
# make sure you are in the vina environment
conda activate vina
scrub.py "smile 1" -o 0-prepared.sdf
mk_prepare_ligand.py -i 0-prepared.sdf -o 0-prepared.pdbqt
```

However, in the folder `Autodock-Vina/ligands`, you will find a file named `ligand-preparation.py` which automatically performs all these iterations for you.

```
# make sure you are in the vina environment
conda activate vina
cd Autodock-Vina/ligands
python ligand-preparation.py
```

This code adds a second column to the file, numbering the molecules from 0 to 300. The output will be saved as `list.csv`. In step 2, the code reads the SMILES strings from the `list.csv` file, adds hydrogens to the molecules, and generates different conformers using molscrub. Then, in step 3, Meeko will create the .pdbqt files for each ligand, saving the files according to their id-num. At the end, you should find the .sdf and .pdbqt files for each SMILES string in your `list.csv` file located in the ligands folder. This is the output of this command:

```
(vina) francesco@Mac ligands % python ligands-preparation.py 
======================================================================
STEP 1: Looking for cheese*.csv file to create list.csv
======================================================================

--- Processing cheese file: ./cheese-full-list-mcule.csv ---
--- Output will be saved to: ./list.csv ---
✓ Created ./list.csv with 'id-num' column added
  Total data rows: 297

======================================================================
STEP 2: Ligand Preparation Pipeline
======================================================================

--- Starting Ligand Preparation Pipeline for: ./list.csv ---

Header row: ['smiles', 'id-num']

[Ligand 1] Processing ID: 0
  SMILES: C1CCC(COc2c3c(nc[nH]3)nc(Nc3ccccc3)n2)CC1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '0-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '0-prepared.pdbqt'

[Ligand 2] Processing ID: 1
  SMILES: COc1ccc(Nc2nc(NCC3CCCO3)c3ccccc3n2)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '1-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '1-prepared.pdbqt'

[Ligand 3] Processing ID: 2
  SMILES: COc1ccc(Nc2nc(NCC3CCCO3)nc3nccnc23)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '2-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '2-prepared.pdbqt'

.......I'M REMOVING SOME OUTPUTS.......

[Ligand 297] Processing ID: 296
  SMILES: CCc1ccc(Nc2nc(NCCCO)nc3c2cnn3C)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '296-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '296-prepared.pdbqt'


======================================================================
--- PIPELINE VERIFICATION REPORT ---
======================================================================

1. FILE EXISTENCE CHECK:
----------------------------------------------------------------------

2. PROCESSING SUMMARY:
----------------------------------------------------------------------
Total ligands processed:          297
Successfully converted to PDBQT:  297
Failed conversions:               0
Missing SDF files:                0
Missing PDBQT files:              0

Success rate:                     100.0%

======================================================================
✓ SUCCESS: All ligands were successfully converted!
======================================================================

--- Pipeline finished: Processed 297 ligands ---
```

# AutoDock Vina

We are now ready to run the Vina docking process. To do this, navigate to the `AutoDock-Vina` folder and execute the following command:

```
# Make sure you are in the vina environment
conda activate vina
python vina-batch.py
```

This command will perform docking for each ligand sequentially, using the following command:

```
vina --receptor receptor/1H1Q-receptor.pdbqt --ligand ligands/0-prepared.pdbqt --config receptor/1H1Q-receptor.box.txt --exhaustiveness 100 --out poses/0-vina-out.pdbqt --num_modes 20
```

Instead of using the Python script `vina-batch.py`, you can also use the `--batch` flag when executing the vina command. However, I prefer to maintain more control over the order of operations and outputs. By using `vina-batch.py`, you will find both the .pdbqt files and the individual docking output files (such as `in-num docking output.txt`) in the `poses` folder. For example, for id-num 0 in `poses/`, we can find the files `0-vina-out.pdbqt` and `0-vina-score.txt`, which look like this:

```
AutoDock Vina v1.2.5
#################################################################
# If you used AutoDock Vina in your work, please cite:          #
#                                                               #
# J. Eberhardt, D. Santos-Martins, A. F. Tillack, and S. Forli  #
# AutoDock Vina 1.2.0: New Docking Methods, Expanded Force      #
# Field, and Python Bindings, J. Chem. Inf. Model. (2021)       #
# DOI 10.1021/acs.jcim.1c00203                                  #
#                                                               #
# O. Trott, A. J. Olson,                                        #
# AutoDock Vina: improving the speed and accuracy of docking    #
# with a new scoring function, efficient optimization and       #
# multithreading, J. Comp. Chem. (2010)                         #
# DOI 10.1002/jcc.21334                                         #
#                                                               #
# Please see https://github.com/ccsb-scripps/AutoDock-Vina for  #
# more information.                                             #
#################################################################

Scoring function : vina
Rigid receptor: receptor/1H1Q-prepared.pdbqt
Ligand: ligands/0-prepared.pdbqt
Grid center: X 6.145 Y 44.175 Z 50.827
Grid size  : X 20 Y 20 Z 20
Grid space : 0.375
Exhaustiveness: 100
CPU: 0
Verbosity: 1

Computing Vina grid ... done.
Performing docking (random seed: 1415709252) ... 
0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
***************************************************

mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -9.723          0          0
   2       -9.464      1.464      4.402
   3       -9.361      2.411      6.166
   4       -9.208     0.5783      1.526
   5       -9.112      1.681       5.27
   6       -9.066        1.5      4.525
   7       -8.966      2.227      5.265
   8       -8.938      2.073      5.325
   9       -8.854      2.638      6.553
  10       -8.823      2.594      6.006
  11       -8.807      2.008      6.655
  12       -8.801       1.46      2.194
  13       -8.748        1.8      4.765
  14       -8.741      2.153      6.883
  15       -8.724      2.381      3.495
  16       -8.722      2.684      5.305
  17       -8.616      1.515        2.4
  18        -8.53      2.398      5.233
  19       -8.495      2.229      6.692
  20       -8.229      2.419      5.247
```

In `poses`, you can also find `ranking.py`, which will take the list of ligands and add one more column with the best docking score from Vina. 

```
python ranking.py
```

We will see a new file called `list_with_affinities.csv` in the folder `Autodock-Vina/poses/`.


### Automatization Vina
To streamline the process, there is a bash script named `run-vina.sh`. Below is the log information generated during the script's execution.
**Note: this is not the actual dataset utilized, but a small sample to illustrate the outputs.**

```
# Make sure you are in the vina environment
conda activate vina
# Make sure that the .sh file is executable
chmod +x run-vina.sh
./run-vina.sh
```

The output will look like this:

```
======================================================================
STEP 1: Looking for cheese*.csv file to create list.csv
======================================================================

--- Processing cheese file: ./cheese-fijijwdk.csv ---
--- Output will be saved to: ./list.csv ---
✓ Created ./list.csv with 'id-num' column added
  Total data rows: 9

======================================================================
STEP 2: Ligand Preparation Pipeline
======================================================================

--- Starting Ligand Preparation Pipeline for: ./list.csv ---

Header row: ['smiles', 'id-num']

[Ligand 1] Processing ID: 0
  SMILES: C1CCC(COc2c3c(nc[nH]3)nc(Nc3ccccc3)n2)CC1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '0-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '0-prepared.pdbqt'

[Ligand 2] Processing ID: 1
  SMILES: COc1ccc(Nc2nc(NCC3CCCO3)c3ccccc3n2)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '1-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '1-prepared.pdbqt'

[Ligand 3] Processing ID: 2
  SMILES: COc1ccc(Nc2nc(NCC3CCCO3)nc3nccnc23)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '2-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '2-prepared.pdbqt'

[Ligand 4] Processing ID: 3
  SMILES: COc1ccc(Nc2nc(NCC3CCCO3)nc3[nH]ncc23)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '3-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '3-prepared.pdbqt'

[Ligand 5] Processing ID: 4
  SMILES: Cc1ccc(Nc2nc(NCC3CCCO3)nc(N)c2[N+](=O)[O-])cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '4-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '4-prepared.pdbqt'

[Ligand 6] Processing ID: 5
  SMILES: Cc1ccc(Nc2nc(NCCCO)c3ccccc3n2)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '5-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '5-prepared.pdbqt'

[Ligand 7] Processing ID: 6
  SMILES: c1ccc(-c2cc(NCC3(Cn4cccn4)CC3)n3nccc3n2)cc1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '6-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '6-prepared.pdbqt'

[Ligand 8] Processing ID: 7
  SMILES: OCCCNc1nc(Nc2ccc(F)cc2)nc2ccccc12
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '7-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '7-prepared.pdbqt'

[Ligand 9] Processing ID: 8
  SMILES: CCCCNc1nc(Nc2ccccc2O)c2cnn(C)c2n1
  > Step 1: Running scrub.py...
  > scrub.py SUCCESS. SDF file created: '8-prepared.sdf'
  > Step 2: Running mk_prepare_ligand.py...
  > mk_prepare_ligand.py SUCCESS. PDBQT file created: '8-prepared.pdbqt'


======================================================================
--- PIPELINE VERIFICATION REPORT ---
======================================================================

1. FILE EXISTENCE CHECK:
----------------------------------------------------------------------

2. PROCESSING SUMMARY:
----------------------------------------------------------------------
Total ligands processed:          9
Successfully converted to PDBQT:  9
Failed conversions:               0
Missing SDF files:                0
Missing PDBQT files:              0

Success rate:                     100.0%

======================================================================
✓ SUCCESS: All ligands were successfully converted!
======================================================================

--- Pipeline finished: Processed 9 ligands ---
Found 9 ligand files
Starting docking process...

Processing ligand 0...
  ✓ Completed: poses/0-vina-score.txt
Processing ligand 1...
  ✓ Completed: poses/1-vina-score.txt
Processing ligand 2...
  ✓ Completed: poses/2-vina-score.txt
Processing ligand 3...
  ✓ Completed: poses/3-vina-score.txt
Processing ligand 4...
  ✓ Completed: poses/4-vina-score.txt
Processing ligand 5...
  ✓ Completed: poses/5-vina-score.txt
Processing ligand 6...
  ✓ Completed: poses/6-vina-score.txt
Processing ligand 7...
  ✓ Completed: poses/7-vina-score.txt
Processing ligand 8...
  ✓ Completed: poses/8-vina-score.txt

======================================================================
DOCKING SUMMARY
======================================================================

Total ligands processed:     9
Successful dockings:         9
Failed dockings:             0

Output files verification:
  Score files (.txt):        9/9
  Pose files (.pdbqt):       9/9

Success rate:                100.0%

All outputs saved in: /Users/francesco/Downloads/terra_quantum/TEST2/Autodock-Vina/poses/
======================================================================

✓ SUCCESS: All ligands docked successfully!
======================================================================

Running from Autodock-Vina directory
Reading ligands data from: ligands/list.csv
Using 'id-num' as the matching column
Found 9 vina-score.txt files
Successfully extracted 9 affinity values

Results saved to: poses/list_with_affinities.csv

============================================================
AFFINITY STATISTICS
============================================================
Total ligands: 9
Ligands with scores: 9
Best affinity: -9.666 kcal/mol
Worst affinity: -8.160 kcal/mol
Mean affinity: -8.822 kcal/mol
Median affinity: -8.826 kcal/mol
============================================================

Top 10 Best Binders:
 id-num  vina_affinity
      0         -9.666
      1         -9.344
      7         -8.856
      5         -8.849
      2         -8.826
      4         -8.623
      6         -8.559
      3         -8.519
      8         -8.160

✓ Processing complete!
```

# Boltz

We can now navigate to the `boltz/` folder, where the `boltz-processing.py` file is located. 

```
# Make sure you are now in boltz-2 environment
conda activate boltz2
python boltz-processing.py              
```

This file can create and run all the other necessary files for us. Alternatively, we can choose to go through each file manually. For each smile, we need to create a .yaml configuration file and then execute it. The `processing.py` script will create .yaml files in the `boltz-configuration-files` folder, where we define the parameters for running the prediction. Specifically, this script will open the file `Autodock-Vina/ligands/list.csv`, which contains the ligands and extracts their SMILES strings. For each SMILES string, it will generate a file named `0.yaml`, `1.yaml`, and so on. In these files, we specify that we will calculate the affinity using the Boltz-2 method, and we have constrained the binding to the pocket where residues 83 and 134 point into the cavity. 

The configuration `0.yaml` looks like this:

```
version: 1  # Optional, defaults to 1
sequences:
  - protein:
      id: A
      sequence: GPLGSMENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREISLLKELNHPNIVKLLDVIHTENKLYLVFEFLHQDLKKFMDASALTGIPLPLIKSYLFQLLQGLAFCHSHRVLHRDLKPQNLLINTEGAIKLADFGLARAFGVPVRTYTHEVVTLWYRAPEILLGCKYYSTAVDIWSLGCIFAEMVTRRALFPGDSEIDQLFRIFRTLGTPDEVVWPGVTSMPDYKPSFPKWARQDFSKVVPPLDEDGRSLLSQMLHYDPNKRISAKAALAHPFFQDVTKPVPHLRL
  - ligand:
      id: B
      smiles: 'C1CCC(COc2c3c(nc[nH]3)nc(Nc3ccccc3)n2)CC1' # Grabbed from the list.csv file, it's every time different
properties:
  - affinity:
      binder: B
constraints:
  - pocket:
      binder: B
      contacts: [ [ A, 83 ], [ A, 134 ] ]
```

After generating this file, the `boltz-processing.py` will run the prediction as follows:

```
boltz predict boltz-configuration-files/0.yaml --use_msa_server --out_dir boltz-results --recycling_steps 1 --sampling_steps 50 --diffusion_samples 3 --step_scale 1.2             
```

For each smile, it will generate a prediction that can be found in `boltz-results/boltz-results_0/predictions/0/affinity.json`, which looks like this:

```
{
    "affinity_pred_value": -0.6497622728347778,
    "affinity_probability_binary": 0.8011776208877563,
    "affinity_pred_value1": -0.6406533718109131,
    "affinity_probability_binary1": 0.7879448533058167,
    "affinity_pred_value2": -0.6588712334632874,
    "affinity_probability_binary2": 0.8144103288650513
}             
```

Once all the ligands are done, you can run the 'boltz-prediction.py' file, 

```
python boltz-prediction.py           
```

which, once again, will take the file with the docking scores from Vina's `list_with_affinities.csv` and it will update it with the Boltz predictions by adding new columns. In the new list, `list_with_affinities_boltz.csv`, we will also find the predictions converted in kcal/mol as explained in their [documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md). So, in the new `list_with_affinities_boltz.csv` there will be the `smiles`, `id-num`, `vina_affinity`, `boltz_affinity_kcalmol`, `avg_affinity_pred_value` (which is the model output score), and finally the `avg_affinity_probability_binary`. 

## Automatization Boltz

Once again, we can avoid running all these steps individually by using the following commands:

```
# Make sure you are now in boltz-2 environment
conda activate boltz2

# Make sure that the .sh file is executable
chmod +x run-boltz.sh
./run-boltz.sh         
```

# Additional Properties

## Sorting

We need to rank the molecules based on their predicted binding energies. Since the scales of the predicted values for Vina and Boltz are similar but with opposite signs, I will calculate an average using the formula `average = (Vina score - Boltz score) / 2`. Alternatively, we could normalize the values and compute the average based on their scales using the formula `average = max model score / (max score - min score)` or any other custom formula. Here I will use the first one, not considering in this sorting the boltz `avg_affinity_probability_binary`.

```
# Make sure you are now in deepchem environment
conda activate deepchem-env
python sorting.py
```

The file `sorting.py` will generate the following output:

```
Total molecules: 297
Molecules with Vina affinity: 297
Molecules with Boltz affinity: 291
Molecules with both affinities: 291

✓ Saved sorted list with 297 molecules to 'list-sorted.csv'
✓ Saved top 10 molecules to 'list-best10.csv'

=== DATA RANGES ===
Vina affinity range: -10.470 to -6.180
Boltz affinity range: 5.254 to 9.476

=== SORTING FORMULA EXPLANATION ===
SIMPLE AVERAGING METHOD:

1. Vina affinity: More negative = better binding (e.g., -9.5 is better than -8.0)
2. Boltz affinity: Higher values = better binding (e.g., 8.0 is better than 6.5)

3. Combined score formula:
   combined_score = (vina_affinity - boltz_affinity_kcalmol) / 2

   Why subtract Boltz?
   - We want to reward LOW vina (e.g., -9.5)
   - We want to reward HIGH boltz (e.g., 8.0)
   - Subtracting high boltz makes the score MORE negative (better)

   Example:
   - Molecule A: Vina=-9.5, Boltz=8.0 → Score=(-9.5-8.0)/2 = -8.75
   - Molecule B: Vina=-8.0, Boltz=6.5 → Score=(-8.0-6.5)/2 = -7.25
   - Molecule A wins (more negative score)

4. Molecules are sorted by combined score (most negative = best)
5. Molecules missing either value are placed at the end

=== TOP 10 MOLECULES ===
                                          smiles  id-num  vina_affinity  boltz_affinity_kcalmol  avg_affinity_pred_value  avg_affinity_probability_binary
0      C1CCC(COc2c3c(nc[nH]3)nc(Nc3ccccc3)n2)CC1       0         -9.723                9.070276                -0.649762                         0.801178
1  c1ccc(CNc2nc(Nc3cc(C4CC4)[nH]n3)c3sccc3n2)cc1      25         -9.314                9.088149                -0.662866                         0.722170
2      COc1ccc(Nc2nc(NCc3cccnc3)nc3[nH]ncc23)cc1     269         -9.008                8.851198                -0.489148                         0.593732
3       Cn1ncc2c(Nc3cc(Cl)ccc3O)nc(NCCCNC=O)nc21      50         -8.790                9.054088                -0.637894                         0.450187
4     NC(=O)c1ccc(Nc2nc(NCC3CCCO3)c3ccccc3n2)cc1      33         -9.701                7.928586                 0.187254                         0.410534
5    Clc1ccc(Nc2nc(NCCC3=CCCCC3)nc3[nH]cnc23)cc1     125         -9.762                7.752512                 0.316340                         0.551613
6        Cc1ccc(Nc2nc(NCCc3ccccc3)nc3nccnc23)cc1     174         -9.970                7.465778                 0.526556                         0.465921
7          CCCCNc1nc(Nc2cc(Cl)ccc2O)c2cnn(C)c2n1      56         -8.543                8.811897                -0.460335                         0.402390
8       CC(C)CCNc1nc(Nc2ccc(F)cc2Cl)c2cnn(C)c2n1     116         -9.512                7.841735                 0.250927                         0.225434
9       CCCCNc1nc(Nc2ccc3nc[nH]c3c2)c2cnn(C)c2n1      27         -8.705                8.642884                -0.336425                         0.268122

=== STATISTICS OF TOP 10 ===
Total molecules with both affinities: 291
Best combined score: -9.397
10th best combined score: -8.674
Average Vina in top 10: -9.303
Average Boltz in top 10: 8.451
Average pred_value in top 10: -0.196
Average probability in top 10: 0.489
```

We will now have a file named `list-sorted.csv`, where all the ligands are ranked based on the average binding energies from the models. Additionally, there is a file called `list-best10.csv` containing the top 10 molecules based on this averaging. 

## RDKit and DeepChem

Finally, we can run the last file, `additional-descriptor.py`

```
python additional-descriptor.py
```

This file is going to create the final table `candidates.csv`. In this table, we will include many descriptors that can be extremely useful for understanding the final candidate molecules. 

Let's begin with the analysis of Lipinski's Rule of Five, which is a fundamental principle in drug discovery used to predict oral bioavailability based on four key molecular properties:

- **Molecular Weight:** Should be ≤500 Da, as larger molecules typically have poor absorption.
- **LogP:** Should be ≤5, since excessive lipophilicity can lead to poor solubility and increased toxicity.
- **Hydrogen Bond Donors:** Should be ≤5, as having too many donors can reduce permeability.
- **Hydrogen Bond Acceptors:** Should be ≤10, as excessive hydrogen bonding can hinder the ability to cross cell membranes. 

This guideline helps in assessing the potential success of a compound as an oral drug. However, we don't know if in this case out molecule is going to be an oral drug. Compounds are considered drug-like if they violate no more than one of these criteria. Additionally, the number of rotatable bonds is assessed, as excessive flexibility (>10 rotatable bonds) can reduce binding affinity due to entropic penalties. So this first section will give us the following columns: `molecular_weight`, `logP`,	`num_h_donors`,	`num_h_acceptors`,	`num_rotatable_bonds`,	`lipinski_violations`,	and `lipinski_pass`.	

The QED (Quantitative Estimation of Drug-likeness) score offers a more detailed evaluation of drug-likeness by combining eight molecular properties, which are weighted according to their significance for oral medications. The score ranges from 0 to 1, with the following classifications:

- 0.7 to 1.0: Excellent drug-likeness
- 0.5 to 0.7: Good drug-likeness
- 0.3 to 0.5: Moderate drug-likeness
- Below 0.3: Poor drug-likeness

These features are added with the following names in the candidates list: `qed_score` and	`qed_classification`.

The Topological Polar Surface Area (TPSA) measures the surface area occupied by polar atoms (oxygen and nitrogen) and their attached hydrogens. This descriptor strongly correlates with membrane permeability and blood-brain barrier penetration. For oral drugs the range should be around 40-100 A^2. Values below 40 A^2 indicate excessive lipophilicity (poor solubility, off-target effects). Values above 100 A^2 suggest poor membrane permeability. 

The ATP binding pocket of kinases contains aromatic residues that engage in π-π stacking interactions with inhibitors.
