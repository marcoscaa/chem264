import numpy as np
from ase.io import read
from ase.units import Ry, Bohr, eV, Angstrom

# Define the magnitude of the electric field in atomic units
efield_magnitude_au = 0.001

# Define the directions of the electric field and corresponding labels
efield_directions = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, 0, 1), (0, 0, -1)]
efield_labels = ['+x', '-x', '+y', '-y', '+z', '-z']

def calculate_born_charges(all_forces, efield_magnitude):
    """
    Calculates Born effective charges from a list of force arrays.

    Args:
        all_forces (list of numpy.ndarray): A list of force arrays, one for each
            electric field direction.  Each array has shape (natoms, 3).
        efield_magnitude (float): The magnitude of the electric field.

    Returns:
        numpy.ndarray: A 3D array of Born effective charges with shape (natoms, 3, 3).
    """
    natoms = all_forces[0].shape[0]
    born_charges = np.zeros((natoms, 3, 3))

    # Calculate Born charges using array slicing
    for atom_index in range(natoms):
        for i in range(3):
            delta_f = (all_forces[2 * i][atom_index] - all_forces[2 * i + 1][atom_index])
            born_charges[atom_index, i, :] = -delta_f / (2 * efield_magnitude)
    return born_charges

# Main script
if __name__ == "__main__":
    # Read atoms objects and get forces
    all_forces = []
    for efield_label in efield_labels:
        output_filename = f'scf_efield_{efield_label}.out'
        try:
            atoms = read(output_filename, format='espresso-out')  # Read the structure
            forces = atoms.get_forces()  # Get the forces in eV/Angstrom
            all_forces.append(forces)
        except Exception as e:
            print(f"Error reading forces from {output_filename}: {e}")
            exit()  # Exit if any file fails to read

    natoms=len(all_forces[0])
    efield_magnitude_V_A = efield_magnitude_au * 36.3609
    born_charges_va = calculate_born_charges(all_forces, efield_magnitude_V_A)
    print("\nBorn Effective Charges (Forces in eV/Angstrom, E-field in Ry/bohr):")
    for atom_index in range(natoms):
        print(f"Atom {atom_index + 1}:")
        for i in range(3):
            for j in range(3):
                print(f"  Z*({i+1},{j+1}) = {-born_charges_va[atom_index, i, j]:.6f}")

    # Write Born effective charges to BORN file
    np.savetxt("BORN", born_charges_va.reshape((-1,9)))
