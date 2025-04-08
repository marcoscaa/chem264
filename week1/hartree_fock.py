import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import time

# --- Configuration Parameters ---
N_points = 401      # Number of grid points (odd number is good for centering)
x_min = -10.0       # Grid minimum (Bohr)
x_max = 10.0        # Grid maximum (Bohr)
R = 1.4             # Internuclear distance (Bohr, equilibrium distance for H2 is ~1.4)
a_soft = 0.01        # Softening parameter for Coulomb potential (prevents division by zero)
convergence_criterion = 1e-7 # Convergence threshold for energy change
max_scf_iterations = 100     # Maximum number of SCF iterations
mixing_alpha = 0.5          # Simple mixing parameter for Fock matrix (0=no mix, 1=full new)

# --- Grid Setup ---
x = np.linspace(x_min, x_max, N_points)
dx = x[1] - x[0]
print(f"Grid setup: {N_points} points from {x_min:.2f} to {x_max:.2f} Bohr (dx = {dx:.4f})")
print(f"Internuclear distance R = {R:.2f} Bohr")
print(f"Coulomb softening parameter a = {a_soft:.2f}")

# --- Define Potentials and Operators ---

# Nuclear positions
R1 = -R / 2.0
R2 = +R / 2.0

# Softened Coulomb interaction function
def soft_coulomb(x1, x2, a):
    """Calculates the softened Coulomb potential 1/sqrt((x1-x2)^2 + a^2)."""
    # Add a small epsilon to denominator inside sqrt for absolute safety, though a^2 should handle it
    # epsilon = 1e-15
    return 1.0 / np.sqrt((x1 - x2)**2 + a**2 ) # + epsilon)

# 1. Nuclear Attraction Potential (V_nuc)
V_nuc = -soft_coulomb(x, R1, a_soft) - soft_coulomb(x, R2, a_soft)
V_nuc_matrix = np.diag(V_nuc) # Diagonal matrix representation

# 2. Kinetic Energy Operator (T) using finite differences
# T = -1/2 * d^2/dx^2
T_matrix = np.zeros((N_points, N_points))
T_matrix = -0.5 * np.diag(np.ones(N_points-1), 1) - 0.5 * np.diag(np.ones(N_points-1), -1) + np.diag(np.ones(N_points), 0)
T_matrix /= (x[1] - x[0])**2  # T in terms of the grid spacing

# Boundary conditions (optional, but often T assumes wavefunction goes to 0 at ends)
# For simplicity, we don't impose strict boundary conditions beyond grid limits here.

# 3. Core Hamiltonian (H_core = T + V_nuc)
H_core = T_matrix + V_nuc_matrix

# 4. Nuclear-Nuclear Repulsion (Constant)
V_nn = soft_coulomb(R1, R2, a_soft) # Use softened version for consistency if R->0
# V_nn = 1.0 / R # Or exact for R > 0

print(f"Nuclear Repulsion V_nn = {V_nn:.6f} Hartree")

# --- SCF Procedure ---

# Initial Guess for the spatial orbital phi
# Let's use a simple normalized Gaussian centered at the origin
# sigma_guess = 1.0
# phi = np.exp(-x**2 / (2 * sigma_guess**2))

# Alternative guess: sum of two Gaussians near nuclei
sigma_guess = 0.8
phi = np.exp(-(x - R1)**2 / (2 * sigma_guess**2)) + \
      np.exp(-(x - R2)**2 / (2 * sigma_guess**2))

# Normalize the initial guess orbital
norm = np.sqrt(np.sum(np.abs(phi)**2) * dx)
phi = phi / norm
print("Initial guess orbital generated and normalized.")

# Previous energy for convergence check
E_total_old = 0.0
E_orbital_old = 0.0
F_old = np.zeros_like(H_core) # For mixing

start_time = time.time()

print("\n--- Starting SCF Iterations ---")
print("Iter | Orbital Energy (Ha) | Total Energy (Ha)  | Delta E (Ha)   | Time (s)")
print("-" * 65)

for iteration in range(max_scf_iterations):
    iter_start_time = time.time()

    # 1. Calculate Electron Density (rho = |phi|^2)
    # For RHF with 2 electrons in orbital phi, density = 2 * |phi|^2
    # But the potential calculation involves integral rho(x') V(x, x'),
    # and for RHF J[phi](x) = integral |phi(x')|^2 V(x, x') dx' phi(x).
    # So we use |phi|^2 to build the potential V_J.
    density = np.abs(phi)**2

    # 2. Calculate Coulomb Potential (V_J) from the density
    V_J = np.zeros_like(x)
    for i in range(N_points):
        # Integral: sum over j of density[j] * interaction(x[i], x[j]) * dx
        V_J[i] = np.sum(density * soft_coulomb(x[i], x, a_soft)) * dx

    V_J_matrix = np.diag(V_J)

    # 3. Construct Fock Operator (F = H_core + V_J)
    # Note: For RHF with 2 electrons in one orbital, J=K, so F = T + V_nuc + J
    F_new = H_core + V_J_matrix

    # Simple Density/Fock matrix mixing (to aid convergence)
    if iteration > 0:
        F = (1 - mixing_alpha) * F_old + mixing_alpha * F_new
    else:
        F = F_new # No mixing on first iteration

    # 4. Solve the Eigenvalue Problem: F * phi = epsilon * phi
    # We need the lowest eigenvalue and corresponding eigenvector
    # scipy.linalg.eigh is good for Hermitian matrices (Fock matrix is)
    # It returns eigenvalues in ascending order
    eigenvalues, eigenvectors = scipy.linalg.eigh(F, b=None) # b=None for standard eigenvalue problem F*v=lambda*v

    # The lowest eigenvalue is the orbital energy
    epsilon = eigenvalues[0]
    # The corresponding eigenvector is the new orbital
    phi_new = eigenvectors[:, 0]

    # 5. Normalize the new orbital
    # Important: Normalization should use the grid spacing dx
    norm = np.sqrt(np.sum(np.abs(phi_new)**2) * dx)
    phi_new = phi_new / norm

    # Ensure phi is real and consistent phase (e.g., positive at first max)
    # Find index of max absolute value
    max_abs_idx = np.argmax(np.abs(phi_new))
    if phi_new[max_abs_idx] < 0:
         phi_new *= -1 # Flip sign if needed

    # Update orbital for next iteration
    phi = phi_new
    F_old = F # Store current Fock matrix for next mixing step

    # 6. Calculate Total Energy
    # E = sum(eps_i) - 1/2 * sum(J_ij - K_ij) + V_nn
    # For RHF with 2 electrons in orbital phi (energy eps):
    # E = 2*eps - Integral[ |phi(x)|^2 * V_J(x) dx ] + V_nn
    E_electron_interaction = np.trapz(density * V_J,x)  # Integral V_J * density dx
    E_total = 2 * epsilon - E_electron_interaction + V_nn

    # 7. Check Convergence
    delta_E = np.abs(E_total - E_total_old)
    iter_end_time = time.time()
    iter_time = iter_end_time - iter_start_time

    print(f"{iteration+1:<4} | {epsilon:<19.8f} | {E_total:<18.8f} | {delta_E:<14.4e} | {iter_time:.2f}")

    if delta_E < convergence_criterion:
        print("-" * 65)
        print(f"SCF converged after {iteration + 1} iterations.")
        break

    # Update old energy for next iteration's comparison
    E_total_old = E_total
    E_orbital_old = epsilon

else: # This else clause executes if the loop completes without break
    print("-" * 65)
    print(f"SCF did not converge within {max_scf_iterations} iterations.")

end_time = time.time()
total_time = end_time - start_time
print(f"Total calculation time: {total_time:.2f} seconds")

# --- Results ---
print("\n--- Final Results ---")
print(f"Converged Orbital Energy (epsilon): {epsilon:.8f} Hartree")
print(f"Converged Total Electronic Energy (E_elec = 2*eps - V_J_int): {E_total - V_nn:.8f} Hartree")
print(f"Converged Total Energy (E_total = E_elec + V_nn): {E_total:.8f} Hartree")

# Calculate final density
final_density = 2 * np.abs(phi)**2 # Physical density (2 electrons)

# --- Plotting ---
plt.figure(figsize=(12, 6))

# Plot Orbital
plt.subplot(1, 2, 1)
plt.plot(x, phi, label=r'$\phi(x)$ (Lowest Energy Orbital)')
plt.title(f'Converged HF Orbital for 1D Hb (R={R:.2f})')
plt.xlabel('Position x (Bohr)')
plt.ylabel('Orbital Amplitude')
plt.axvline(R1, color='r', linestyle='--', linewidth=0.8, label=f'Nucleus 1 ({R1:.2f})')
plt.axvline(R2, color='r', linestyle='--', linewidth=0.8, label=f'Nucleus 2 ({R2:.2f})')
plt.grid(True, alpha=0.5)
plt.legend()
plt.tight_layout()

# Plot Density
plt.subplot(1, 2, 2)
plt.plot(x, final_density, label=r'$\rho(x) = 2|\phi(x)|^2$')
plt.title(f'Electron Density for 1D Hb (R={R:.2f})')
plt.xlabel('Position x (Bohr)')
plt.ylabel('Electron Density')
plt.axvline(R1, color='r', linestyle='--', linewidth=0.8, label=f'Nucleus 1 ({R1:.2f})')
plt.axvline(R2, color='r', linestyle='--', linewidth=0.8, label=f'Nucleus 2 ({R2:.2f})')
plt.grid(True, alpha=0.5)
plt.legend()
plt.tight_layout()

plt.suptitle(f"1D Hartree-Fock Hb (Grid Method) - E_total = {E_total:.6f} Ha")
plt.subplots_adjust(top=0.9) # Adjust title position
plt.show()
