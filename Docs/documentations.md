# Documentation

## Table of Contents
* [Detailed Documentation](#detailed-documentation)
* [Files Preparation](#files-preparation)
  * [Receptor Preparation](#receptor-preparation)
  * [Ligand Preparation](#ligands-preparation)
* [AutoDock Vina](#autodock-vina)
* [Boltz](#boltz)

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


### run-vina.sh automatization
To streamline the process, `Autodock-Vina` includes a bash script named `run-vina.sh`. Below is the log information generated during the script's execution.
**Note: this is not the actual dataset utilized, but a small sample to illustrate the outputs.**

```
# Make sure you are in the vina environment
(base) francesco@Mac % conda activate vina
(vina) francesco@Mac % ./run-vina.sh 
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

Now you can navigate in boltz folder where again you can find the `boltz-processing.py` file that can run all the files for you. Otherwise, you can manually go over each of them. 

The `processing.py` will place in the `boltz-configuration-files` folder the .yaml files, where we define how we want to run the prediction. Specifically, this file will open `Autodock-Vina/ligands/list.csv`, which contains the ligands, and extract the SMILES strings. For each SMILES string, it will generate a file named `0.yaml`, `1.yaml`, and so on. In these files, we specify that we will calculate the affinity using the Boltz-2 method, and we have constrained the binding within the pocket where residues 83 and 134 point into the cavity. 

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

After generating this file, the `boltz-processing.py` will run the prediction as following:

```
$ boltz predict boltz-configuration-files/0.yaml --use_msa_server --out_dir boltz-results --recycling_steps 1 --sampling_steps 50 --diffusion_samples 3 --step_scale 1.2
                
```

For each smile, it will generate a prediction that can be found in boltz-results/boltz-results_0/predictions/0/affinity.json, which looks like this:

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

Once all the ligands are done, you can run the 'boltz-processing.py' file, which, once again, will take the file with the docking scores from Vina (list_with_affinities.csv) and it will update it with the Boltz prediction by adding another column. In the new list, list_with_affinities_boltz.csv, you will find the prediction already converted in kcal/mol as explained in their [documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).
