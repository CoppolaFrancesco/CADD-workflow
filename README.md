# CADD-workflow
A fully automated workflow for Drug Discovery

# Introduction
This workflow automates the process of screening multiple molecules against a specific protein target. Once the user selects the dataset of molecules and the receptor, two independent docking protocols will be executed: Vina docking and Boltz. A dedicated protocol will merge the two docking results, and additional features will be added. In the end, the most suitable candidates are presented to the user. 

All the outputs are provided for a dataset of ≈300 molecules for an alternative inhibitor of CDK2. There are two ways to utilize this workflow: one involves automated bash.sh commands for both Vina and Boltz, along with extracting several candidates with additional properties. The second would be to follow the documentation file and run each step separately to have precise control over each step. In both cases, the commands (.py or .sh) must be executed with a previously activated environment, as outlined below.  

# Environments and Installations

I recommend separate environments for the 3 steps: [AutoDock Vina](https://autodock-vina.readthedocs.io/en/latest/installation.html#python-bindings-linux-and-mac-only), [Boltz-2](https://github.com/forlilab/molscrub), and the additional features. The reason is that they might have conflicts in the same environment, and this also allows us to have more control in each step of this workflow. In the future, we may have a single environment and a more harmonized workflow structure. 

## Autodock-Vina

```
conda create -n vina python=3
conda activate vina
```

The following command to install NumPy and AutoDock Vina:

```
conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
pip install vina
```

In the same environment of vina we also need to install [molscrub](https://github.com/forlilab/molscrub).

```
git clone git@github.com:forlilab/molscrub.git
cd molscrub
pip install -e .
```

We also need [Meeko](https://meeko.readthedocs.io/en/release-doc/) for generating the .pdbqt files.

```
conda install meeko
```

## For Boltz-2: 

We also need a dedicated environment for boltz 

```
conda create -n boltz2 python=3.11 -y
conda activate boltz2

pip install boltz
```

In my case, I don't have a dedicated GPU, and I will run everything on the CPU.

## Deepchem: additional properties

I will also need a dedicated environment for the DeepChem package used for calculating additional properties:

```
conda create -n deepchem_env
conda activate deepchem_env
pip install deepchem
```


# Documentation
# Receptor preparation

For this demonstration, I'm using the CDK2 (PDB 1H1Q). Inside these files, there are usually many waters (some of them are important!), cofactors, other proteins, and ligands. For now, for simplicity, we are going to take only the receptor. I recommend two ways: using [VMD](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) or the PDB Reader & Manipulator in [CHARMM-GUI](https://www.charmm-gui.org/?doc=input/pdbreader). In step 1, you will be prompted to select which chain you want to model; in my case, it’s chain A (segid PROA). Step 2 is particularly important because it allows us to manipulate the PDB file. We need to carefully consider the protonation states and any missing residues. Once done, we can download the pdb file, which I also provide in the folder data/receptor/1H1Q-CHARMM-GUI.pdb along with the following generated files.

We need to generate the PDBQT file for our receptor, which is a modified .pdb file that includes partial charges and atom types. 

In my case, I'm using [Meeko](https://meeko.readthedocs.io/en/release-doc/):
```
$ conda install meeko
```

After installation, we can generate the PDBQT file of our receptor, along with the docking pocket box, which I approximated based on the crystallographic structure of the complex with the ligand NU6094.
```
mk_prepare_receptor.py -i 1H1Q-CHARMM-GUI.pdb -o 1H1Q-prepared -p -v --box_size 20 20 20 --box_center 6.145 44.175 50.827 
```
This command should now generate the following files 1H1Q-prepared.pdbqt, 1H1Q-prepared.box.txt and 1H1Q-prepared.box.pdb which concluded our receptor preparation.

# Ligand preparation

For the dataset preparation, I used the [cheese search](https://cheese.deepmedchem.com/), in particular by looking only in the Mcule Full database, by using a search method called [Shape: 3D volume overlap](https://chemrxiv.org/engage/chemrxiv/article-details/67250915f9980725cfcd1f6f). In the dataset I tried to have few analogs of the ligand (NU6094) and some diverse compounds in both structural similarity and number of
rotatable bonds. I selected around 300 molecules and exported as .csv file. The .csv will have now the following structure: first column is the smile, second column the id-num 0 to 300. Now in the folder data/ligands you will also find a file called ligand-preparation.py which if you run it 

```
$ cd data/ligands
$ python ligand-preparation.py
```

Will first (step 1) read the smiles from the .csv file add the hydrogens to the molecules and also generate different conformers with the help of molscrub, and then in step 2 with Meeko is going to generate the .pdbqt files for each ligands saving the files by their resIDs. At the end, in the ligands folder you should see the .sdf and .pdbqt files for each smiles in your list.csv.

# AutoDock Vina docking

Now we are ready to run the Vina docking. To do so in the folder AutoDock-Vina/ just run the following command:

```
$ python vina-batch.py
```

This command will run for each ligand the docking one by one using this command:

```
$ vina --receptor receptor/1H1Q-receptor.pdbqt --ligand ligands/0-prepared.pdbqt --config receptor/1H1Q-receptor.box.txt --exhaustiveness 100 --out poses/0-vina-out.pdbqt --num_modes 20
```

Alternatively to the python vina-batch.py one can also use the flag --batch when running the vina command, but I prefer to have more control over the order and the outputs. In fact, by using my command in the poses/ folder, you can find both the pdbqt files and the specific docking output. In poses, you can also find ranking.py, which will take the list of ligands and add one more column with the best docking score from Vina. This new file is called list_with_affinities.csv and can be found with all the rest of the docking .pdbqt files in the folder Autodock-Vina/poses/.

Here is the 

```
(vina) francesco@Mac TEST2 % ./run-vina.sh 
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

Header row: ['smiles', 'id-num', 'id-num', 'id', 'database', 'db_id', 'similarity', 'caco2_wang', 'clearance_hepatocyte_az', 'clearance_microsome_az', 'half_life_obach', 'ld50_zhu', 'lipophilicity_astrazeneca', 'ppbr_az', 'solubility_aqsoldb', 'vdss_lombardo', 'ames', 'bbb_martins', 'bioavailability_ma', 'cyp2c9_substrate_carbonmangels', 'cyp2c9_veith', 'cyp2d6_substrate_carbonmangels', 'cyp2d6_veith', 'cyp3a4_substrate_carbonmangels', 'cyp3a4_veith', 'dili', 'herg', 'hia_hou', 'pgp_broccatelli', 'molecular_weight', 'formal_charge', 'clogp', 'heavy_atoms', 'h_bond_acceptors', 'h_bond_donor', 'rotatable_bonds', 'num_of_rings', 'molar_refractivity', 'number_of_atoms', 'topological_surface_area_mapping']

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
(vina) francesco@Mac TEST2 % 
```

# ML docking: Boltz

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
