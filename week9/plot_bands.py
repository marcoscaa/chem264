import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(3.3,3.3))

# load data
data = np.loadtxt('ge_bands.dat.gnu')

#High-symmetry points (see 02_pp_bands.out) 
hs_pts = {0:'$\Gamma$',1.0:'X',1.3536:'U',1.9659:'L',2.8320:'$\Gamma$',3.8926:'K',4.2462:'W',4.7462:'X'}

# Define K-points (k) and eigenvalues (bands)
k = np.unique(data[:, 0])
bands = np.reshape(data[:, 1], (-1, len(k)))

# Plot each band (\omega(k))
for band in range(len(bands)):
    plt.plot(k, bands[band, :], linewidth=1, alpha=0.5, color='b')
plt.xlim(min(k), max(k))
plt.ylim(bottom=-5)

# Fermi energy
plt.axhline(8.7649,  linewidth=0.75, color='k', alpha=0.5)

# High symmetry k-points 
for K in list(hs_pts.keys()):
    plt.axvline(K, linewidth=0.75, color='k', alpha=0.25)

# text labels
plt.xticks(ticks=list(hs_pts.keys()), labels=list(hs_pts.values()))
plt.ylabel("Energy (eV)")
plt.text(0.8, 7.8, 'Fermi energy', fontsize='small')
plt.tight_layout()
plt.show()
