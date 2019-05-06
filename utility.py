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

MAYACHEMTOOLS_DIR = os.path.abspath("mayachemtools")

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
        
        with open(self.sdf_filepath + ".csv", 'r') as f:
            for line in f.readlines():
                if "Cmpd" in line:
                    line = line.split(';')[5].replace('"','')
                    features = [int(i) for i in line.split(" ")]

        # Clean up the temporary files
        shutil.rmtree(temp_dir)
        return features


if __name__=="__main__":
    # Example: Extracting TPATF features
     ft = FeatureGenerator()
     ft.load_smiles("O=C(Cc1ccccc1)Nc2ncc(s2)C3CCC3")
     features = ft.extract_tpatf()
