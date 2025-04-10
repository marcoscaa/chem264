from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from ase.io import write
from sys import argv

def generate_pwscf_input_from_cif(structure, output_filename="pwscf.in"):

    reader = CifParser(structure)
    structure: Structure = reader.get_structures(primitive=True)[0]  # Get the first structure

    # Convert Pymatgen Structure to ASE Atoms object
    atoms: Atoms = structure.to_ase_atoms()

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
            'mixing_mode': 'plain',
        },
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


if __name__ == "__main__":
    cif_filename = argv[1] 
    output_file = "scf.in"  # You can change the output filename

    generate_pwscf_input_from_cif(cif_filename, output_file)

