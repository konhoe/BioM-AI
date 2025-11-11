#!/usr/bin/env python3
"""
CIF to Slab PDB Converter with Supercell Expansion
Converts a CIF file to a slab structure with supercell expansion and exports as PDB format
"""

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
import argparse
import numpy as np


def create_slab_with_supercell(cif_file, miller_index=(0, 0, 1), min_slab_size=10.0, 
                                min_vacuum_size=15.0, supercell=(1, 1, 1),
                                output_pdb="slab_structure.pdb"):
    """
    Create a slab structure from a CIF file with supercell expansion and save as PDB
    
    Parameters:
    -----------
    cif_file : str
        Path to the input CIF file
    miller_index : tuple
        Miller indices for the slab surface (default: (0,0,1))
    min_slab_size : float
        Minimum slab thickness in Angstroms (default: 10.0)
    min_vacuum_size : float
        Minimum vacuum layer thickness in Angstroms (default: 15.0)
    supercell : tuple
        Supercell expansion in (a, b, c) directions (default: (1,1,1))
    output_pdb : str
        Output PDB filename (default: "slab_structure.pdb")
    """
    
    # Read the CIF file
    print(f"Reading CIF file: {cif_file}")
    parser = CifParser(cif_file)
    structure = parser.parse_structures()[0]
    
    print(f"\nOriginal structure:")
    print(f"  Formula: {structure.composition.reduced_formula}")
    print(f"  Lattice parameters: a={structure.lattice.a:.3f}, "
          f"b={structure.lattice.b:.3f}, c={structure.lattice.c:.3f} Å")
    print(f"  Number of atoms: {len(structure)}")
    
    # Generate slab
    print(f"\nGenerating slab with Miller index {miller_index}")
    print(f"  Min slab size: {min_slab_size} Å")
    print(f"  Min vacuum size: {min_vacuum_size} Å")
    
    slabgen = SlabGenerator(
        structure,
        miller_index,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        center_slab=True,
        primitive=False
    )
    
    # Get all possible slabs and select the first one
    slabs = slabgen.get_slabs()
    
    if not slabs:
        raise ValueError("No slabs generated. Try different parameters.")
    
    print(f"\nFound {len(slabs)} possible slab terminations")
    slab = slabs[0]
    
    print(f"\nSlab structure (before supercell):")
    print(f"  Formula: {slab.composition.reduced_formula}")
    print(f"  Lattice parameters: a={slab.lattice.a:.3f}, "
          f"b={slab.lattice.b:.3f}, c={slab.lattice.c:.3f} Å")
    print(f"  Number of atoms: {len(slab)}")
    
    # Create supercell
    if supercell != (1, 1, 1):
        print(f"\nCreating supercell with expansion {supercell}")
        slab.make_supercell(supercell)
        print(f"  New lattice parameters: a={slab.lattice.a:.3f}, "
              f"b={slab.lattice.b:.3f}, c={slab.lattice.c:.3f} Å")
        print(f"  Number of atoms: {len(slab)}")
    
    # Convert to ASE Atoms object
    adaptor = AseAtomsAdaptor()
    atoms = adaptor.get_atoms(slab)
    
    # Write custom PDB format with HETATM records
    write_custom_pdb(atoms, output_pdb)
    
    print(f"\nSlab structure saved to: {output_pdb}")
    
    return slab


def write_custom_pdb(atoms, filename):
    """
    Write atoms to PDB file with HETATM format and proper formatting
    
    Parameters:
    -----------
    atoms : ase.Atoms
        ASE Atoms object
    filename : str
        Output filename
    """
    
    with open(filename, 'w') as f:
        # Write CRYST1 record
        cell = atoms.get_cell()
        a, b, c = cell.lengths()
        alpha, beta, gamma = cell.angles()
        
        f.write(f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
                f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1\n")
        
        f.write("MODEL     1\n")
        
        # Group atoms by element
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # Sort by element and then by z-coordinate
        element_order = {}
        for i, sym in enumerate(symbols):
            if sym not in element_order:
                element_order[sym] = []
            element_order[sym].append((i, positions[i]))
        
        # Write HETATM records
        atom_index = 1
        for element in sorted(element_order.keys()):
            atom_list = element_order[element]
            # Sort by z, then y, then x coordinate
            atom_list.sort(key=lambda x: (x[1][2], x[1][1], x[1][0]))
            
            for orig_idx, pos in atom_list:
                x, y, z = pos
                
                # Format: HETATM with proper spacing
                # Residue name is element + "4" (like TI4)
                res_name = f"{element:>2s}4"
                
                line = (f"HETATM{atom_index:5d} {element:>2s}   {res_name:3s} Z"
                       f"{atom_index:4d}    "
                       f"{x:8.3f}{y:8.3f}{z:8.3f}"
                       f"  1.00 -0.00          {element:>2s}  \n")
                
                f.write(line)
                atom_index += 1
        
        f.write("ENDMDL\n")
    
    print(f"  Written {atom_index-1} atoms in HETATM format")


def main():
    parser = argparse.ArgumentParser(
        description="Convert CIF file to slab structure in PDB format with supercell expansion"
    )
    parser.add_argument(
        "cif_file",
        help="Input CIF file path"
    )
    parser.add_argument(
        "-m", "--miller",
        nargs=3,
        type=int,
        default=[0, 0, 1],
        help="Miller indices for slab surface (default: 0 0 1)"
    )
    parser.add_argument(
        "-s", "--slab-size",
        type=float,
        default=10.0,
        help="Minimum slab thickness in Angstroms (default: 10.0)"
    )
    parser.add_argument(
        "-v", "--vacuum-size",
        type=float,
        default=15.0,
        help="Minimum vacuum thickness in Angstroms (default: 15.0)"
    )
    parser.add_argument(
        "-x", "--supercell",
        nargs=3,
        type=int,
        default=[1, 1, 1],
        help="Supercell expansion (a b c), e.g., 10 8 1 (default: 1 1 1)"
    )
    parser.add_argument(
        "-o", "--output",
        default="slab_structure.pdb",
        help="Output PDB filename (default: slab_structure.pdb)"
    )
    
    args = parser.parse_args()
    
    miller_index = tuple(args.miller)
    supercell = tuple(args.supercell)
    
    try:
        create_slab_with_supercell(
            args.cif_file,
            miller_index=miller_index,
            min_slab_size=args.slab_size,
            min_vacuum_size=args.vacuum_size,
            supercell=supercell,
            output_pdb=args.output
        )
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())