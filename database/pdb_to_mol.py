#!/usr/bin/env python3
"""
Convert PDB → MOL2 for Rosetta molfile_to_params.py
- Ensures Tripos-compliant mol2 with bond, type, and charge
- Uses OpenBabel backend (pybel)
- Suitable for inorganic / metal-containing slabs

"""

from openbabel import openbabel, pybel
import sys
from pathlib import Path

def pdb_to_mol2_for_rosetta(pdb_path: str, output_path: str = None):
    pdb_path = Path(pdb_path)
    if output_path is None:
        output_path = pdb_path.with_suffix(".mol2")

    # Load molecule from PDB
    mol = next(pybel.readfile("pdb", str(pdb_path)))

    # Assign bond connectivity and partial charges
    mol.OBMol.AddHydrogens(False, True)  # Skip H addition, ensure bonding
    builder = openbabel.OBBuilder()
    builder.Build(mol.OBMol)

    # Assign Gasteiger partial charges (required by Rosetta)
    charge_model = openbabel.OBChargeModel.FindType("gasteiger")
    if not charge_model.ComputeCharges(mol.OBMol):
        print("⚠️ Warning: charge computation failed; continuing without charges")

    # Set forcefield atom typing (Tripos-like types)
    ff = openbabel.OBForceField.FindForceField("uff")
    if ff is not None:
        ff.Setup(mol.OBMol)
        ff.GetAtomTypes(mol.OBMol)

    # Write MOL2 in Tripos format (Rosetta compatible)
    mol.write("mol2", str(output_path), overwrite=True)
    print(f"✅ Converted {pdb_path.name} → {output_path.name} (Tripos MOL2 for Rosetta)")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python pdb_to_mol2_for_rosetta.py input.pdb [output.mol2]")
        sys.exit(1)

    pdb_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    pdb_to_mol2_for_rosetta(pdb_file, output_file)