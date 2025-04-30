import numpy as np
from ase.io import read
from sys import argv

#Unit conversion
eVtoJ=1.60218e-19
A2toM2=1.E-20

#Read QE output file
bulk = read(argv[1], index=-1)
surface = read(argv[2], index=-1)

#Energy per atom of bulk
Ebulk = bulk.get_potential_energy()
Ebulk /= len(bulk.get_atomic_numbers())

#Energy of the slab
Eslab = surface.get_potential_energy()
N = len(surface.get_atomic_numbers())

#Surface area
cell=surface.get_cell()
area=np.linalg.norm(np.cross(cell[0],cell[1]))

#Surface energy
Esurf = Eslab - N*Ebulk
Esurf/= 2*area
Esurf *= eVtoJ/A2toM2

print(Esurf)

