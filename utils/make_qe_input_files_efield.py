from ase.io import read, write
import numpy as np
from sys import argv

# Define the magnitude of the electric field
efield_magnitude_au = 0.001

# Define the directions of the electric field
efield_directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
efield_labels = ['+x', '-x', '+y', '-y', '+z', '-z']

# Name of the Quantum Espresso output file from the vc-relax calculation
output_file = argv[1]

try:
    # Read the structure from the output file
    atoms = read(output_file, format='espresso-out')
except FileNotFoundError:
    print(f"Error: Output file '{output_file}' not found.")
    exit()
except Exception as e:
    print(f"Error reading output file: {e}")
    exit()

# Create input files for SCF calculations with external electric fields
for i, direction in enumerate(efield_directions):
    input_filename = f'scf_efield_{efield_labels[i]}.in'

    # Basic input dictionary (modify as needed for your system)
    input_data = {
        'control': {
            'calculation': 'scf',
            'restart_mode': 'from_scratch',
            'prefix': 'scf_efield',
            'pseudo_dir': './',  # Adjust if your pseudo directory is different
            'outdir': './out',        # Adjust if you want a different output directory
            'lelfield': True,      # Enable external electric field
            'disk_io': 'none',
            'tprnfor': True,
        },
        'system': {
            'ibrav': 0,
            'nat': len(atoms),
            'ntyp': len(set(atoms.get_chemical_symbols())),
            'ecutwfc': 60.0,      # Adjust cutoff energy as needed
        },
        'electrons': {
            'conv_thr': 1.0e-8,   # Adjust convergence threshold as needed
            'mixing_beta': 0.7,   # Adjust mixing parameter as needed
            'diagonalization': 'david',
            'efield_cart(1)': efield_magnitude_au * direction[0],
            'efield_cart(2)': efield_magnitude_au * direction[1],
            'efield_cart(3)': efield_magnitude_au * direction[2],
            'startingpot': 'file',
        },
        'ions': {
            'ion_dynamics': 'none', # No ionic relaxation for SCF
        },
        'cell': {
            'cell_dynamics': 'none', # No cell relaxation for SCF
        },
    }

    # Add atomic positions and cell from the relaxed structure
    input_data['atomic_positions'] = atoms.get_chemical_symbols()
    input_data['atomic_positions_type'] = 'crystal'  # Assuming crystal coordinates
    input_data['cell_parameters'] = atoms.cell.cellpar()
    pseudopot = {"O":"O.upf","Ge":"Ge.upf"}
    kpts=(6,4,4)

    # Write the input file
    try:
        write(input_filename, atoms, format='espresso-in', input_data=input_data, pseudopotentials=pseudopot, kpts=kpts)
        print(f"Generated input file: {input_filename} with E-field along {efield_labels[i]} direction.")
    except Exception as e:
        print(f"Error writing input file '{input_filename}': {e}")

print("Finished generating SCF input files with external electric fields.")
