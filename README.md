## Dependencies
```bash
pip install rdkit
pip install py3dmol
```

## Tips
The functions take SMILES to make the following type of figures:

1. Simple Pubchem style plot saved as png with clear background.
2. Colorful stick and ball or stick-only 3D views, not able to save. I use screenshot and then use Powerpoint to remove background.

How to get SMILES string? 

1. From xyz file, you can use rdkit.
2. Go to Pubchem website and input the molecule name or formula.
3. Ask ChatGPT but usually not accurate.
4. Go to the Pubchem Sketch website and draw the molecule, then copy the SMILES.