# CADD-workflow
A fully automated workflow for Drug Discovery

# Introduction
I will write it later but VINA + BOLTZ2

# Environments

I recommend separate environments, one for [AutoDock Vina ](https://autodock-vina.readthedocs.io/en/latest/installation.html#python-bindings-linux-and-mac-only) and one for [Boltz-2](https://github.com/forlilab/molscrub). 

For Vina just run in the terminal:

```
$ conda create -n vina python=3
$ conda activate vina
```

The following command to install NumPy and AutoDock Vina:

```
$ conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme
$ pip install vina
```

For Bolz-2: 

```
$ conda create boltz_conda
$ conda activate boltz_conda
```

Now you can just download and install it. In my case, I don't have a dedicated GPU and I will run everything on the CPU 

```
git clone https://github.com/jwohlwend/boltz.git
cd boltz
pip install -e .
```

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

Before run this in the Vina enviroment we need to install [molscrub](https://github.com/forlilab/molscrub).

```
git clone git@github.com:forlilab/molscrub.git
cd molscrub
pip install -e .
```

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
