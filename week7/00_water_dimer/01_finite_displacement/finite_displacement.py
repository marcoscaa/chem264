import os
from ase.io import read, write
from sys import argv
import numpy as np

def generate_finite_differences(qe_output_file, d=0.01):
    """
    Generates Quantum ESPRESSO input files for finite difference calculations
    based on a relaxed structure from a Quantum ESPRESSO output file.

    Args:
        qe_output_file (str): Path to the Quantum ESPRESSO output file
            from a relaxation calculation (e.g., 'relax.out').
        d (float, optional): The displacement in Angstroms. Defaults to 0.01.
    """
    try:
        # 1. Read the last configuration from the Quantum ESPRESSO output file.
        #    ASE automatically takes the last structure for a geometry optimization.
        atoms = read(qe_output_file, format='espresso-out')
    except Exception as e:
        print(f"Error reading the Quantum ESPRESSO output file: {e}")
        return

    # Get the original cell
    cell = atoms.get_cell()

    # 2. Loop over each atom and each Cartesian coordinate (x, y, z).
    for i, atom in enumerate(atoms):
        for j in range(3):
            # Create a copy of the original atomic positions.
            positions = atoms.get_positions().copy()

            # Apply the displacement to the jth coordinate of the ith atom.
            positions[i, j] += d

            # Create a new Atoms object with the displaced positions.
            displaced_atoms = atoms.copy()  # Start with a copy to avoid modifying original
            displaced_atoms.set_positions(positions)
            displaced_atoms.set_cell(cell) # keep the same cell

            # 3. Generate a Quantum ESPRESSO input file for the displaced structure.
            input_filename = f"atom_{i+1}_dir_{j+1}.in"
            #  Use a dictionary to define the input parameters.  This is cleaner
            #  and less error-prone than string concatenation.
            input_params = {
                'control': {
                    'calculation': 'scf',  # Perform an SCF calculation
                    'restart_mode': 'from_scratch',
                    'prefix': f'atom_{i+1}_dir_{j+1}', # unique prefix 
                    'pseudo_dir': './',  
                    'outdir': f'./tmp_{i+1}_{j+1}', # unique outdir
                    'tprnfor': True,
                    'tstress': False,
                    'disk_io': 'none',
                },
                'system': {
                    'ecutwfc': 60.0,  
                    'nat': len(atoms),
                    'ntyp': len(set(atoms.get_chemical_symbols())), # number of unique atom types
                },
                'electrons': {
                    'conv_thr': 1.0e-8,  
                    'mixing_mode': 'plain',
                    'mixing_beta': 0.5,
                },
            }

            pseudo = {"O":"O.upf","H":"H.upf"}

            write(input_filename, displaced_atoms, format='espresso-in',
                  kpts=None,  # Example k-points, adjust as needed
                  pseudopotentials=pseudo,
                  input_data=input_params)

    print("Finished generating Quantum ESPRESSO input files.")

if __name__ == "__main__":
    qe_output_file = argv[1]
    generate_finite_differences(qe_output_file, d=0.01)

