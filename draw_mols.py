from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

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


def draw_3d_stick_ball(smiles_string, stick_scale=0.1, sphere_scale=0.25):
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
    p.setStyle({'stick': {'radius': stick_scale}, 'sphere': {'scale': sphere_scale}})
    p.zoomTo()
    p.show()
    return p

def draw_molecule_in_black(smiles, save_name="molecule.png"):
    # Create the molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    rdDepictor.Compute2DCoords(mol)  # Generate 2D coordinates
    
    # Set up drawing options for black and white
    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)  # Create a Cairo image
    options = drawer.drawOptions()
    options.useBWAtomPalette()  # Set the palette to black and white
    options.setBackgroundColour((1,1,1,0))
    # Draw the molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    # Save and display the image
    img = drawer.GetDrawingText()
    with open(save_name, "wb") as f:
        f.write(img)