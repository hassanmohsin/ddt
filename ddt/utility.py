# Utility script for feature generation
# Md Mahmudulla Hassan
# The University of Texas at El Paso
# Last Modified: 05/13/2019
from __future__ import print_function, absolute_import
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import shutil
import subprocess
import errno
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions
import base64
import numpy as np

current_dir = os.path.dirname(os.path.realpath(__file__))
MAYACHEMTOOLS_DIR = os.path.join(current_dir, "mayachemtools")


class SmilesToImage:
    def __init__(self, smiles):
        self.smiles = smiles
        self.temp_dir = tempfile.mkdtemp()
        self.png_file = os.path.join(self.temp_dir, "mol.png")

    def toPNG(self, get_binary=False, output=None):
        if output is not None:
            self.png_file = output
        # Set the drawing options
        # DrawingOptions.atomLabelFontSize = 55
        # DrawingOptions.dotsPerAngstrom = 100
        # DrawingOptions.bondLineWidth = 3.0

        # Convert the SMILES into a mol object
        m = Chem.MolFromSmiles(self.smiles)
        # Calculate the coordinates
        AllChem.Compute2DCoords(m)
        # Draw the mol
        img = Draw.MolToImage(m, size=(300, 300))
        img.save(self.png_file)
        # Draw.MolToFile(m, self.svg_file)
        # Convert the svg to png (for high quality image)
        # cairosvg.svg2png(url=self.svg_file, write_to=self.png_file)
        # Convert into binary and return
        if get_binary:
            binary_image = None
            with open(self.png_file, "rb") as f:
                binary_image = base64.b64encode(f.read())
                shutil.rmtree(self.temp_dir)

            return binary_image


class SDFToImage:
    def __init__(self, sdf):
        self.sdf = sdf
        self.temp_dir = tempfile.mkdtemp()
        self.png_file = os.path.join(self.temp_dir, "mol.png")

    def toPNG(self, get_binary=False, output=None, max=8, molsPerRow=4):
        if output is not None:
            self.png_file = output
        # Set the drawing options
        # DrawingOptions.atomLabelFontSize = 55
        # DrawingOptions.dotsPerAngstrom = 100
        # DrawingOptions.bondLineWidth = 3.0

        # Convert the SDF into a mol object
        suppl = Chem.SDMolSupplier(self.sdf)
        ms = [x for x in suppl if x is not None]
        for m in ms:
            tmp = AllChem.Compute2DCoords(m)
        # Draw the mol
        img = Draw.MolsToGridImage(
            ms[:max],
            molsPerRow=molsPerRow,
            subImgSize=(200, 200),
            legends=[x.GetProp("_Name") for x in ms[:max]],
        )
        img.save(self.png_file)
        # Convert into binary and return
        if get_binary:
            binary_image = None
            with open(self.png_file, "rb") as f:
                binary_image = base64.b64encode(f.read())
                shutil.rmtree(self.temp_dir)

            return binary_image


class FeatureGenerator:
    def __init__(self):
        self.sdf_filepath = None
        self.smiles = None

    def smiles_to_sdf(self, smiles):
        # Try to get the rdkit mol
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("ERROR: Invalid SMILES")
            return None
        # Compute 2D coordinates
        AllChem.Compute2DCoords(mol)
        mol.SetProp("smiles", smiles)
        temp_dir = tempfile.mkdtemp()
        sdf_filepath = os.path.join(temp_dir, "temp.sdf")
        w = Chem.SDWriter(sdf_filepath)
        w.write(mol)
        w.flush()
        return sdf_filepath

    def load_smiles(self, smiles):
        self.smiles = smiles
        self.sdf_filepath = self.smiles_to_sdf(smiles)

    def load_sdf(self, sdf_filepath):
        self.sdf_filepath = sdf_filepath

    def load_smi_file(
        self, smiles_file, titleLine, smilesColumn, nameColumn, delimiter
    ):
        """
        smiles_file (str): file containing SMILES strings
        smilesColumn (int): index of the smiles column
        titleLine (bool): True if there is a header in the file, False if not
        delimiter (str): delimiter of the file, e.g., ',' for csv
        """
        temp_dir = tempfile.mkdtemp()
        self.sdf_filepath = os.path.join(temp_dir, "temp.sdf")
        mols = AllChem.SmilesMolSupplier(
            smiles_file,
            delimiter=delimiter,
            smilesColumn=smilesColumn,
            nameColumn=nameColumn,
            titleLine=titleLine,
        )
        w = Chem.SDWriter(self.sdf_filepath)
        counter = 0
        for mol in mols:
            if not mol:
                continue
            AllChem.Compute2DCoords(mol)
            w.write(mol)
            counter = counter + 1
        w.flush()

        print(f"{len(mols)}/{counter} compounds are converted into sdf")

    def extract_tpatf(self, csv_file=None):
        features = []
        script_path = os.path.join(
            MAYACHEMTOOLS_DIR, "bin/TopologicalPharmacophoreAtomTripletsFingerprints.pl"
        )
        # Generate the TPATF features
        # Check if the sdf file exists
        if not os.path.isfile(self.sdf_filepath):
            print("SDF file not found")
            return None

        # Extract tpatf features
        temp_dir = tempfile.mkdtemp()
        temp_file = os.path.join(temp_dir, "temp")
        # Use below command to use the compound name as compound id (selleck compounds sdf has it)
        # command = "perl " + script_path + " -r " + temp_file + " --DataFieldsMode Specify --DataFields Name --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + self.sdf_filepath
        # Use below command to use the available compound id in the sdf file (most sdf files have it)
        command = (
            "perl "
            + script_path
            + " -r "
            + temp_file
            + " --CompoundIDMode MolnameOrLabelPrefix --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o "
            + self.sdf_filepath
        )
        # command = "perl " + script_path + " -r " + temp_file + " --DataFieldsMode CompoundID --CompoundIDMode MolnameOrLabelPrefix --CompoundID Cmpd --CompoundIDLabel MolID --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + self.sdf_filepath
        os.system(command)
        output_csv = temp_file + ".csv"
        if not os.path.isfile(output_csv):
            print("ERROR: TPATF features wasn't extracted")
            return None

        if csv_file:
            shutil.copy(output_csv, os.path.join(os.getcwd(), csv_file))

        compound_list = []
        content = []

        # TODO: Add warning with line number that fails to load because of invalid decode error
        with open(output_csv, "rb") as f:
            for line in f:
                content.append(line.decode("utf-8", "ignore"))

        content = [c.replace('"', "") for c in content]  # remove (") from the content
        content = [
            c.split(";") for c in content
        ]  # split features from the other information

        # Separate the compound features
        for con in content[1:]:  # First item doesn't have features
            compound_list.append(con[0].split(",")[0])
            features.append([int(i) for i in con[-1].split(" ")])

        # with open(output_csv, 'r') as f:
        #    for line in f.readlines():
        #        #if "Cmpd" in line:
        #        if line[1:5] == "Cmpd":
        #            #compound_list.append(list[1:
        #            line = line.split(';')[5].replace('"','')
        #            features.append([int(i) for i in line.split(" ")])

        # Clean up the temporary files
        shutil.rmtree(temp_dir)
        return compound_list, np.array(features).reshape((-1, 2692)).astype(np.float32)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Feature generation from SMILES and SDF"
    )
    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument(
        "--smiles", action="store", dest="smiles", help='SMILES string as "SMILES"',
    )
    input.add_argument("--sdf", action="store", dest="sdf", help="SDF file location")
    input.add_argument(
        "--smiles_file", action="store", dest="smiles_file", help="SMILES file"
    )
    parse_file = parser.add_argument_group("Argument for SMILES file parsing")
    parse_file.add_argument(
        "--title",
        action="store_true",
        dest="title",
        help="If the file contains a header",
    )
    parse_file.add_argument(
        "--smiles_column",
        action="store",
        dest="smiles_column",
        help="Column index for SMILES entries",
    )
    parse_file.add_argument(
        "--name_column",
        action="store",
        dest="name_column",
        help="Column index for SMILES names",
    )
    parse_file.add_argument(
        "--delimiter",
        action="store",
        dest="delimiter",
        help="Delimiter, e.g., ' ' or ','",
    )

    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument(
        "--npy", action="store", dest="npy", help="Numpy file for writing the features"
    )
    output.add_argument(
        "--csv", action="store", dest="csv", help="CSV file to write the features",
    )
    args = parser.parse_args()

    ft = FeatureGenerator()

    if args.smiles:
        ft.load_smiles(args.smiles)
    elif args.smiles_file:
        ft.load_smi_file(
            args.smiles_file,
            titleLine=args.title,
            smilesColumn=int(args.smiles_column),
            delimiter=args.delimiter,
            nameColumn=int(args.name_column),
        )
    elif args.sdf:
        ft.load_sdf(args.sdf)
    else:
        raise NotImplementedError

    if args.csv:
        _, _ = ft.extract_tpatf(csv_file=args.csv)
    else:
        compounds, features = ft.extract_tpatf()
        np.save(args.npy, features)
        with open(args.npy + ".txt", "w") as f:
            f.writelines(",".join([c for c in compounds]))
