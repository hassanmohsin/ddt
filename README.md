# ddt
A toolkit for drug discovery research

## Usage:

### Terminal
```
 $ python -m ddt.utility --help
usage: utility.py [-h]
                  (--smiles SMILES | --sdf SDF | --smiles_file SMILES_FILE)
                  (--npy | --csv)

Feature generation from SMILES and SDF

optional arguments:
  -h, --help            show this help message and exit
  --smiles SMILES       SMILES string as "SMILES"
  --sdf SDF             SDF file location
  --smiles_file SMILES_FILE
                        SMILES file
  --npy                 Save features to .npy
  --csv                 Get the features in a csv file
  ```

  For a smiles file, prepare the file as below:
  
  ```
  NAME,SMILES
  levobupivacaine,CCCCN1CCCC[C@H]1C(=O)NC1=C(C)C=CC=C1C
  (S)-nicardipine,COC(=O)C1=C(C)NC(C)=C([C@H]1C1=CC(=CC=C1)[N+]([O-])=O)C(=O)OCCN(C)CC1=CC=CC=C1
  (S)-nitrendipine,CCOC(=O)C1=C(C)NC(C)=C([C@@H]1C1=CC(=CC=C1)[N+]([O-])=O)C(=O)OC
  levdobutamine,C[C@@H](CCC1=CC=C(O)C=C1)NCCC1=CC=C(O)C(O)=C1
  aminopterin,NC1=NC2=NC=C(CNC3=CC=C(C=C3)C(=O)N[C@@H](CCC(O)=O)C(O)=O)N=C2C(N)=N1
  aminomethylbenzoic acid,NCC1=CC=C(C=C1)C(O)=O
  phenylbutanoic acid,OC(=O)CCCC1=CC=CC=C1
  azacitidine,NC1=NC(=O)N(C=N1)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O
  fluorouracil,FC1=CNC(=O)NC1=O
  ```

  This is a CSV file. First and seconds column contain the compound names and the SMILES respectively. If the compound names are not available, use dummy texts.

### Script

```
Python 3.6.10 |Anaconda, Inc.| (default, Jan  7 2020, 21:14:29)
[GCC 7.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from ddt.utility import FeatureGenerator
>>> ft = FeatureGenerator()
>>> ft.load_smiles("CCCC")
>>> c, f = ft.extract_tpatf()

TopologicalPharmacophoreAtomTripletsFingerprints.pl: Starting...

Processing options...
Checking input SD file(s)...

Processing file /tmp/tmpwvy0hglv/temp.sdf...
Generating text file /tmp/tmp2mc2qv61/temp.csv...

Number of compounds: 1
Number of compounds processed successfully during fingerprints generation: 1
Number of compounds ignored during fingerprints generation: 0

TopologicalPharmacophoreAtomTripletsFingerprints.pl:Done...

Total time:  0 wallclock secs ( 0.05 usr +  0.00 sys =  0.05 CPU)

>>> print(c)
['Cmpd1']
>>> f.shape
(1, 2692)

```