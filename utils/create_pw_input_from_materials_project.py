from ase import Atoms
from ase.io import write
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Structure
from sys import argv

def generate_pwscf_input_from_mp(structure, output_filename="pwscf.in"):

    atoms = AseAtomsAdaptor.get_atoms(structure)

    input_params = {
        'control': {
            'calculation':'scf',
            'restart_mode': 'from_scratch',
            'prefix': 'ge',
            'pseudo_dir': './',
            'outdir': './out',
            'tprnfor': True,
            'tstress': True,
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
        }
    }

    pseudopotentials={}
    for symbol in atoms.get_chemical_symbols():
        if symbol not in pseudopotentials:
            pseudopotentials[symbol] = symbol+".upf"
    kpoints_grid_config = "888"

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
                return structure
            else:
                print(f"No structure found for MP-ID: {mp_id}")
                return None
    except Exception as e:
        print(f"Incorrect MP-ID '{mp_id}'")
        return None


if __name__ == "__main__":
    # Replace with the desired Materials Project ID and your API key
    material_id_to_extract = argv[1] 
    your_api_key = "qQccZYnGqnS0ad2qoZ92chnoQdBvKgSb"  # Replace with your actual API key
    output_file = "scf.in"  # You can change the output filename

    #Authenticate your Materials Project ID
    structure = get_structure_from_mp(your_api_key,material_id_to_extract)

    generate_pwscf_input_from_mp(structure, output_file)

