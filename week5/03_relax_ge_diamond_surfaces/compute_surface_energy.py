import numpy as np
from ase.io import read
from sys import argv

#Read QE output file
bulk = read(argv[1], index=-1)
surface = read(argv[2], index=-1)

#Energy per atom of bulk
Ebulk = bulk.get_potential_energy()
Ebulk /= len(bulk.get_atomic_numbers())

#Energy of the slab
Eslab = surface.get_potential_energy()
N = len(surface.get_atomic_numbers())

#Surface energy
Esurf = Eslab - N*Ebulk
Esurf/= 2

print(Esurf)

