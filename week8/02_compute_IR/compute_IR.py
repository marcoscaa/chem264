import numpy as np
from ase.units import _c
from phonopy import Phonopy, load

e_charge = 1.602176634e-19  # C
amu_to_kg = 1.66053906660e-27  # kg/amu
Bohr=5.29177e-11 #m/Bohr

#Load the Phonopy input file
phonon = load("phonopy_disp.yaml",force_sets_filename='FORCE_SETS')

# Extract the Born effective charges
born_charges = np.loadtxt('BORN').reshape((-1,3,3))

print("Calculating IR intensities...")
frequencies,eigenvectors = phonon.get_frequencies_with_eigenvectors((0,0,0))
primitive_cell = phonon.primitive
volume = np.linalg.det(primitive_cell.primitive_matrix) * Bohr**3
natoms = len(primitive_cell)


# Convert Born charges to proper units (C)
born_charges_e_A = born_charges * e_charge

ir_intensities = [] 
# Loop over q-points (which should be Gamma only for IR)
# In Phonopy, the Gamma point is the first q-point (index 0)
q = 0  # Fix q to be the Gamma point
intensity = 0.0
supercell_atoms = phonon.supercell # Get the atoms of the supercell
print(len(supercell_atoms))
for j in range(len(frequencies)):  # Loop over phonon modes
    # Calculate the polarization vector
    e_j = eigenvectors[:, j]  # Eigenvector for mode j at Gamma
    polarization_vector = np.zeros(3, dtype=complex)
    for a in range(natoms):  # Iterate over atoms in primitive cell
        for d in range(3):
            polarization_vector += (born_charges_e_A[a, d, :] * e_j[a])

    # Sum the squares of the polarization vector components
    intensity = np.sum(np.abs(polarization_vector)**2)
    # Convert intensity to Debye^2/amu
    ir_intensities.append(intensity * 2*np.pi * Bohr**2 * natoms / volume / _c / amu_to_kg / 3.)

# Print the IR spectrum to the console.
print("----------------------------------------------------------------")
print(" Frequency (THz)   Intensity (SI)   Wavenumber (cm-1)")
print("----------------------------------------------------------------")

frequencies_factor=1
if ir_intensities is not None:
    for freq, intensity in zip(frequencies, ir_intensities):
        if np.iscomplex(freq):
            freq = 0.0
        wavenumber = freq * 29.9792458  # Convert THz to cm-1
        print(f"{freq * frequencies_factor:10.4f}   {intensity:14.4e}   {wavenumber:12.4f}")
else:
    for freq in frequencies:
         if np.iscomplex(freq):
            freq = 0.0
         wavenumber = freq * 29.9792458
         print(f"{freq * frequencies_factor:10.4f}   {0.0:14.4f}   {wavenumber:12.4f}")
print("----------------------------------------------------------------")

# Write the IR spectrum to a file, if requested.
print(f"Writing IR spectrum to file: IR.dat")
with open("IR.dat", 'w') as f:
    f.write("# Frequency (THz)  Intensity (SI)  Wavenumber (cm-1)\n")
    for freq, intensity in zip(frequencies, ir_intensities):
        if np.iscomplex(freq):
            freq = 0.0
        wavenumber = freq * 29.9792458
        f.write(f"{freq * frequencies_factor:10.4f}  {intensity:14.4e}  {wavenumber:12.4f}\n")
