# Utility script for feature generation
# Md Mahmudulla Hassan
# The University of Texas at El Paso
# Last Modified: 05/05/2019

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
import cairosvg

current_dir = os.path.dirname(os.path.realpath(__file__))
MAYACHEMTOOLS_DIR = os.path.join(current_dir, "mayachemtools")

def smiles_to_sdf(smiles):
    # Try to get the rdkit mol
    mol = Chem.MolFromSmiles(smiles)
    # Compute 2D coordinates
    AllChem.Compute2DCoords(mol)
    mol.SetProp("smiles", smiles)
    temp_dir = tempfile.mkdtemp()
    sdf_filepath = os.path.join(temp_dir, "temp.sdf")
    w = Chem.SDWriter(sdf_filepath)
    w.write(mol)
    w.flush()
    print(sdf_filepath)
    return sdf_filepath

class SmilesToImage:
    def __init__(self, smiles):
        self.smiles = smiles
        self.temp_dir = tempfile.mkdtemp()
        self.png_file = os.path.join(self.temp_dir, "mol.png")
        self.svg_file = os.path.join(self.temp_dir, "mol.svg")

    def toPNG(self):
        # Set the drawing options
        DrawingOptions.atomLabelFontSize = 55
        DrawingOptions.dotsPerAngstrom = 100
        DrawingOptions.bondLineWidth = 3.0
        
        # Conver the SMILES into a mol object
        m = Chem.MolFromSmiles(self.smiles)
        # Calculate the coordinates
        AllChem.Compute2DCoords(m)
        # Draw the mol
        Draw.MolToFile(m, self.svg_file)
        # Convert the svg to png (for high quality image)
        cairosvg.svg2png(url=self.svg_file, write_to=self.png_file)
        # Convert into binary and return
        binary_image = None
        with open(self.png_file, "rb") as f:
            binary_image = base64.b64encode(f.read())
            shutil.rmtree(self.temp_dir)

        return binary_image


class FeatureGenerator:
    
    def __init__(self):
        self.sdf_filepath = None
        self.smiles = None
        
    def load_smiles(self, smiles):
        self.smiles = smiles
        self.sdf_filepath = smiles_to_sdf(self.smiles)
    
    def load_sdf(self, sdf_filepath):
        self.sdf_filepath = sdf_filepath    
    
    def extract_tpatf(self):
        features = []
        script_path = os.path.join(MAYACHEMTOOLS_DIR, "bin/TopologicalPharmacophoreAtomTripletsFingerprints.pl")
        # Generate the TPATF features
        # Check if the sdf file exists
        if not os.path.isfile(self.sdf_filepath):
            print("SDF file not found")
            return None
        
        # Extract tpatf features
        temp_dir = tempfile.mkdtemp()    
        temp_file = os.path.join(temp_dir, "temp")
        command = "perl " + script_path + " -r " + temp_file + " --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + self.sdf_filepath
        os.system(command)

        with open(temp_file + ".csv", 'r') as f:
            for line in f.readlines():
                #if "Cmpd" in line:
                if line[1:5] == "Cmpd":
                    line = line.split(';')[5].replace('"','')
                    features.append([int(i) for i in line.split(" ")])

        # Clean up the temporary files
        shutil.rmtree(temp_dir)
        return features


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="An utility script for small tasks.")
    parser.add_argument('--smiles', action='store', dest='smiles', required=False, help='SMILES string as "SMILES"')
    parser.add_argument('--sdf', action='store', dest='sdf', required=False, help='SDF file location')
    parse_dict = vars(parser.parse_args())
    
    # Example: Extracting TPATF features
    ft = FeatureGenerator()
    if parse_dict['smiles'] is not None:
        ft.load_smiles(parse_dict['smiles'])
    elif parse_dict['sdf'] is not None:
        ft.load_sdf(parse_dict['sdf'])
        features = ft.extract_tpatf()
        print(len(features))
