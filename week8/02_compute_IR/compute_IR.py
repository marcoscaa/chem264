import numpy as np
from ase.units import Debye, Angstrom, electron_mass, V_Angstromsq
from phonopy import Phonopy, load
from phonopy.file_IO import write_vasp_poscar, write_FORCE_CONSTANTS
from phonopy.interface.calculator import get_default_physical_units

#Load the Phonopy input file
phonon = load("phonopy_disp.yaml")

# Get the forces from the 'FORCES' files.
print("Reading forces from calculations...")
force_sets = phonon.get_force_sets()

# Extract the Born effective charges
born_charges = np.loadtxt('BORN').reshape((-1,3,3))

# Calculate the force constants.
print("Calculating force constants...")
phonon.set_force_sets(force_sets)
phonon.produce_force_constants()

# Calculate the phonon frequencies and eigenvectors.
print("Calculating phonon frequencies and eigenvectors...")
phonon.calculate_frequencies()

print("Calculating IR intensities...")
frequencies = phonon.get_frequencies()
eigenvectors = phonon.get_eigenvectors()
primitive_cell = phonon.get_primitive_cells()
volume = primitive_cell.get_volume()

# Convert Born charges to proper units (e * Angstrom)
born_charges_e_A = born_charges / electron_mass * (Angstrom * electron_mass)

ir_intensities = [] 
# Loop over q-points (which should be Gamma only for IR)
# In Phonopy, the Gamma point is the first q-point (index 0)
q = 0  # Fix q to be the Gamma point
intensity = 0.0
for j in range(len(frequencies[q])):  # Loop over phonon modes
    # Calculate the polarization vector
    e_j = eigenvectors[q, :, j]  # Eigenvector for mode j at Gamma
    polarization_vector = np.zeros(3, dtype=complex)
    for a, atom in enumerate(primitive_cell):
        for d in range(3):
            polarization_vector[d] += (born_charges_e_A[a, d, :] * e_j[primitive_cell.get_global_atom_index(a)])

    # Sum the squares of the polarization vector components
    intensity += np.sum(np.abs(polarization_vector)**2)
# Convert intensity to Debye^2/amu
ir_intensities.append(intensity * V_Angstromsq / volume)

# Print the IR spectrum to the console.
print("\nIR Spectrum (Frequencies in THz, Intensities in Debye^2*amu^-1):")
print("----------------------------------------------------------------")
print(" Frequency (THz)   Intensity (Debye^2/amu)   Wavenumber (cm-1)")
print("----------------------------------------------------------------")

frequencies = phonon.get_frequencies()
if ir_intensities is not None:
    for freq, intensity in zip(frequencies, ir_intensities):
        if np.iscomplex(freq):
            freq = 0.0
        wavenumber = freq * 29.9792458  # Convert THz to cm-1
        print(f"{freq * frequencies_factor:10.4f}   {intensity:14.4f}   {wavenumber:12.4f}")
else:
    for freq in frequencies:
         if np.iscomplex(freq):
            freq = 0.0
         wavenumber = freq * 29.9792458
         print(f"{freq * frequencies_factor:10.4f}   {0.0:14.4f}   {wavenumber:12.4f}")
print("----------------------------------------------------------------")

# Write the IR spectrum to a file, if requested.
if output_filename:
    print(f"Writing IR spectrum to file: {output_filename}")
    with open(output_filename, 'w') as f:
        f.write("# Frequency (THz)  Intensity (Debye^2/amu)  Wavenumber (cm-1)\n")
        for freq, intensity in zip(frequencies, ir_intensities):
            if np.iscomplex(freq):
                freq = 0.0
            wavenumber = freq * 29.9792458
            f.write(f"{freq * frequencies_factor:10.4f}  {intensity:14.4f}  {wavenumber:12.4f}\n")
    print("Done.")
