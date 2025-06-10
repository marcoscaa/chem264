from ase.io import read, write
from glob import glob
from sys import argv

# outcar contanining one (single-point) or multiple (trajectory) DFT frames
in_filename = argv[1]
out_filename = 'input.extxyz'

# read all frames into a list of ase.Atoms objects
all_atoms = read(in_filename, format='espresso-out', index=":")
write(out_filename, all_atoms, format='extxyz')
