from ase.io import write, read
from ase.build import find_optimal_cell_shape, make_supercell
from sys import argv

def orthogonalize(non_ortho_atoms):
    #Compute transformation matrix P
    P = find_optimal_cell_shape(non_ortho_atoms.get_cell(),target_size=2,target_shape='sc')
    #Make supercell - orthogonal
    ortho_atoms = make_supercell(non_ortho_atoms,P)
    #Return orthogonal cell atoms object
    return ortho_atoms

def generate_pwscf_input(atoms, output_filename="pwscf.in"):

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
    non_orthogonal_atoms = read(argv[1], format='espresso-out', index=-1) 
    orthogonal_atoms = orthogonalize(non_orthogonal_atoms)
    output_file = "scf.in"  # You can change the output filename

    generate_pwscf_input(orthogonal_atoms, output_file)

