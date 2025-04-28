from ase import Atoms
from ase.io import write
from ase.build import surface, make_supercell
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
from sys import argv

def generate_pwscf_input_from_mp(atoms, output_filename="pwscf.in"):

    input_params = {
        'control': {
            'calculation':'relax',
            'restart_mode': 'from_scratch',
            'prefix': 'ge',
            'pseudo_dir': './',
            'outdir': './out',
            'tprnfor': True,
            'tstress': False,
        },
        'system': {
            'ibrav': 0,  # 0 for general cell
            'ecutwfc':  60.0,
            'occupations': 'smearing',
            'smearing': 'gaussian',
            'degauss': 0.01,
        },
        'electrons': {
            'mixing_beta': 0.4,
        },
        'ions': {
            'ion_dynamics': 'bfgs',
        }
    }

    pseudopotentials={}
    for symbol in atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            pseudopotentials[symbol] = symbol+".upf"
    kpoints_grid_config = "441"

    write(output_filename, 
            atoms, 
            format='espresso-in', 
            pseudopotentials=pseudopotentials,
            kpts=tuple(kpoints_grid_config),
            input_data=input_params)

def get_structure_from_mp(api_key, mp_id):
    try:
        with MPRester(api_key) as mpr:
            structure = mpr.get_structure_by_material_id(mp_id)
            if structure:
                analyzer = SpacegroupAnalyzer(structure)
                conventional = analyzer.get_conventional_standard_structure()
                return conventional
            else:
                print(f"No structure found for MP-ID: {mp_id}")
                return None
    except Exception as e:
        print(f"Incorrect MP-ID '{mp_id}'")
        return None

def make_surface(structure, miller_index):

    #Convert PyMatGen structure to ase atoms object
    atoms = AseAtomsAdaptor.get_atoms(structure)

    #Create a surface with miller index miller_index
    #The slab has 4 layers, and 10 A vacuum between periodic replicas
    slab = surface(atoms,miller_index,layers=2,vacuum=5,periodic=True)

    return slab

def string_to_tuple_map(s):
  """Converts a 3-digit string to a tuple of integers using the map() function."""
  if not isinstance(s, str) or len(s) != 3 or not s.isdigit():
    return "Invalid input: String must be 3 digits."
  return tuple(map(int, s))

def randomly_displace(atoms, amplitude=0.1):
    """
    Randomly displaces the atoms in an ASE Atoms object.

    Args:
        atoms (ase.Atoms): The Atoms object to displace.
        amplitude (float): The maximum magnitude of the displacement in Angstroms.
                           Default is 0.1 Angstroms.

    Returns:
        ase.Atoms: A new Atoms object with the displaced atomic positions.
    """
    new_atoms = atoms.copy()  # It's good practice to work on a copy
    n_atoms = len(new_atoms)
    displacements = np.random.uniform(low=-amplitude, high=amplitude, size=(n_atoms, 3))
    new_atoms.positions += displacements
    return new_atoms

if __name__ == "__main__":
    # Replace with the desired Materials Project ID and your API key
    material_id_to_extract = argv[1] 
    miller_index = string_to_tuple_map(argv[2]) #should be a 3-digits integer, ex: 100
    na=int(argv[3]) #Replicate x direction na times
    nb=int(argv[4]) #Replicate y direction nb times
    your_api_key = "qQccZYnGqnS0ad2qoZ92chnoQdBvKgSb"  # Replace with your actual API key
    output_file = "surface_relax.in"  # You can change the output filename

    #Authenticate your Materials Project ID
    structure = get_structure_from_mp(your_api_key,material_id_to_extract)

    #Make surface slab with miller index miller_index
    atoms_slab = make_surface(structure, miller_index)

    #Make surface supercell
    atoms_slab = make_supercell(atoms_slab, [[na,0,0],[0,nb,0],[0,0,1]])

    #Randomly displace atoms to move than slightly away from perfect crystal
    atoms_slab = randomly_displace(atoms_slab,0.02) 

    generate_pwscf_input_from_mp(atoms_slab, output_file)

