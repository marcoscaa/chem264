import matplotlib.pyplot as plt
import numpy as np

# load data
energy, dos, idos = np.loadtxt('ge_dos.dat', unpack=True)

#Define fermi energy (see 00_scf.out)
fermi_e = 8.7649 

# make plot
plt.figure(figsize = (3, 3))
plt.plot(energy, dos, linewidth=0.75, color='blue')
plt.yticks([])
plt.xlabel('Energy (eV)')
plt.ylabel('DOS')
plt.axvline(x=fermi_e, linewidth=0.5, color='k')
plt.ylim(bottom=0)
plt.fill_between(energy, 0, dos, where=(energy < fermi_e), facecolor='blue', alpha=0.25)
plt.text(fermi_e-1, 1.2, 'Fermi energy', fontsize='small', rotation=90)
plt.xlim([-5,20])
plt.tight_layout()
plt.show()
