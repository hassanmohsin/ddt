# ddt
A toolkit for drug discovery research

## Usage:

### Terminal
```
 $ python -m ddt.utility --help
usage: utility.py [-h]
                  (--smiles SMILES | --sdf SDF | --smiles_file SMILES_FILE)
                  [--title] [--smiles_column SMILES_COLUMN]
                  [--name_column NAME_COLUMN] [--delimiter DELIMITER]
                  (--npy NPY | --csv CSV)

Feature generation from SMILES and SDF

optional arguments:
  -h, --help            show this help message and exit
  --smiles SMILES       SMILES string as "SMILES"
  --sdf SDF             SDF file location
  --smiles_file SMILES_FILE
                        SMILES file
  --npy NPY             Numpy file for writing the features
  --csv CSV             CSV file to write the features

Argument for SMILES file parsing:
  --title               If the file contains a header
  --smiles_column SMILES_COLUMN
                        Column index for SMILES entries
  --name_column NAME_COLUMN
                        Column index for SMILES names
  --delimiter DELIMITER
                        Delimiter, e.g., ' ' or ','
  ```

Example: `python -m ddt.utility --smiles_file oo.csv --title --smiles_column 3 --name_column 1 --delimiter , --npy test`

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