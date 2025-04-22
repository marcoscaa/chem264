from ase import Atoms
from ase.build import bulk
import numpy as np

def generate_germanium_diamond_xyz(filename="germanium_diamond.xyz"):
    """
    Generates the FCC diamond crystal structure of Germanium using ASE and saves it in XYZ format.

    Args:
        filename (str, optional): The name of the file to save the structure to.
            Defaults to "germanium_diamond.xyz".
    """
    # Create the FCC diamond structure for Germanium.
    # The 'diamond' structure in ASE is the same as FCC with a two-atom basis.
    germanium_diamond = bulk("Ge", crystalstructure="diamond", a=5.657) # a=5.657 is the lattice constant in Angstrom

    # Center the atoms in the cell.  This is often helpful for visualization.
    germanium_diamond.center()

    # Save the structure in XYZ format.
    germanium_diamond.write(filename, format="extxyz")

    print(f"Germanium diamond structure saved to {filename}")

def generate_germanium_conventional_xyz(filename="germanium_conventional.xyz"):
    """
    Generates the conventional FCC cell of Germanium using ASE and saves it in XYZ format.

    Args:
        filename (str, optional): The name of the file to save the structure to.
            Defaults to "germanium_conventional.xyz".
    """
    # Lattice constant
    a = 5.657  # Angstrom

    # Create the conventional FCC cell directly
    conventional_cell = Atoms('Ge8',
                               cell=[[a, 0, 0], [0, a, 0], [0, 0, a]],  # Cell vectors for FCC
                               scaled_positions=[
                                   (0.00, 0.00, 0.00),
                                   (0.50, 0.50, 0.00),
                                   (0.50, 0.00, 0.50),
                                   (0.00, 0.50, 0.50),
                                   (0.25, 0.25, 0.25),  # Positions for the 4 atoms in the basis
                                   (0.75, 0.75, 0.25),
                                   (0.75, 0.25, 0.75),
                                   (0.25, 0.75, 0.75)
                               ],
                               pbc=True)

    # Center the atoms in the cell.
    conventional_cell.center()

    # Save the structure in XYZ format.
    conventional_cell.write(filename, format="extxyz")

    print(f"Germanium conventional cell structure saved to {filename}")

if __name__ == "__main__":
#    generate_germanium_diamond_xyz()
     generate_germanium_conventional_xyz()
