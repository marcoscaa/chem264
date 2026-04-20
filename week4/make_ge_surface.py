from ase import Atoms
from ase.build import bulk, surface
import numpy as np

def generate_germanium_slab(
    layers=3,
    vacuum=10.0,
    index=(0,0,1),
    filename="germanium_slab.xyz",
):
    """
    Generates a Ge(001) slab model using ASE and saves it in XYZ format.

    Args:
        layers (int, optional): The number of atomic layers in the slab. Defaults to 3.
        vacuum (float, optional): The amount of vacuum (in Angstroms) to add
            on both sides of the slab. Defaults to 10.0.
        filename (str, optional): The name of the file to save
            the structure to. Defaults to "germanium_001_slab.xyz".
    """
    # Lattice constant of Germanium
    a = 5.657  # Angstrom

    # Create the bulk FCC diamond structure
    bulk_germanium = bulk("Ge", crystalstructure="diamond", a=a, cubic=True)

    # Create the (001) surface slab
    # The (0, 0, 1) miller index corresponds to the z-direction
    slab = surface(
        bulk_germanium,
        indices=index,  # Miller indices
        layers=layers,
        vacuum=vacuum / 2.0, # Vacuum on one side. ASE adds on both sides.
    )

    slab = delete_undercoordinated_atoms(slab,"Ge","Ge",2.6,2)

    # Center the atoms in the cell.
    slab.center()

    # Save the structure in XYZ format
    slab.write(filename, format="extxyz")

def delete_undercoordinated_atoms(atoms,atype,atype2,rcut,mincoord):
    ind_atoms = [ atom.index for atom in atoms if atom.symbol == atype2 ]
    coord_number = np.full(len(atoms.get_chemical_symbols()),10)
    for ind, atom in enumerate(atoms):
        if atom.symbol == atype:
            distances = atoms.get_distances(atom.index,ind_atoms,mic=True)
            filt=np.logical_and(distances<rcut,distances>0.1)
            coord_number[ind] = len(distances[filt])

    del atoms[coord_number<mincoord]
    return atoms

if __name__ == "__main__":
    # Generate a 3-layer slab with 10 Å vacuum
    generate_germanium_slab(layers=5, vacuum=10.0, index=(1,1,1))
