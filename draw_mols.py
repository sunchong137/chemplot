from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol

def draw_3d_sticks(smiles_string):
    # Convert SMILES to RDKit Molecule object
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol)  # Add hydrogen atoms
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry with UFF force field
    
    # Convert RDKit molecule to a 3D structure viewable by py3Dmol
    mb = Chem.MolToMolBlock(mol)
    p = py3Dmol.view(width=400, height=400)
    p.addModel(mb, 'sdf')
    p.setStyle({'stick': {}})
    p.zoomTo()
    p.show()
    return p


def save_molecule_png(smiles_string, file_name="molecule_3D_view.png"):
    """
    Function to generate and save a 3D plot of a molecule from a SMILES string to a PNG file.
    
    Parameters:
    - smiles_string (str): The SMILES string representing the molecule.
    - file_name (str): The output file name for the PNG image (default: molecule_3D_view.png).
    """
    
    # Generate RDKit molecule object from SMILES
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol)  # Add hydrogen atoms
    
    # Generate 3D coordinates and optimize geometry
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)  # Optimize geometry with UFF force field
    
    # Convert RDKit molecule to a 3D structure viewable by py3Dmol
    mb = Chem.MolToMolBlock(mol)
    p = py3Dmol.view(width=400, height=400)
    p.addModel(mb, 'sdf')
    p.setStyle({'stick': {}})
    p.zoomTo()
    p.show()
    # Save the 3D view as a PNG file
    png_data = p.png()  # Capture the PNG binary data
    with open(file_name, "wb") as f:
        f.write(png_data)

    print(f"Molecule saved to {file_name}")