"""
mean_field.py — Hartree (mean-field) SCF for 1D H2.

Approximation
-------------
Each electron moves in the mean electrostatic field of the TOTAL electron
density rho(x) = 2|phi(x)|^2.  This includes each electron's interaction
with itself (self-interaction error).

The Hartree potential is

    V_H[rho](x) = integral rho(x') / sqrt((x-x')^2 + a^2) dx'

and each electron satisfies the single-particle equation

    [T + V_en + V_H[rho]] phi = eps phi

The density is updated self-consistently using DIIS (Direct Inversion of the
Iterative Subspace) acceleration on the density residual:

    r_n = rho_n - rho_from_diag(F[rho_n])

DIIS finds the optimal linear combination of past densities that minimises
the residual norm, which breaks the oscillatory 2-cycle that arises from the
strong self-interaction in the Hartree equations.

Total energy with the correct double-counting correction (factor 1/2):

    E = 2*eps - (1/2) * integral rho(x) V_H[rho](x) dx + V_nn

Comparison with Hartree-Fock (see hartree_fock.py):
  The Hartree energy E = 2*h_11 + 2*J_11 + V_nn includes an extra J_11 term
  (self-interaction error).  In HF the exchange operator exactly cancels this,
  giving E_HF = 2*h_11 + J_11 + V_nn < E_Hartree.

Usage
-----
    python mean_field.py
"""

import numpy as np
import scipy.linalg as la

from potentials import v_en, v_nn, v_hartree
from utils import make_grid, kinetic_matrix, normalize
from plot import plot_orbital_and_density

# ── Parameters ────────────────────────────────────────────────────────────────
R        = 1.4    # Internuclear distance (a.u.)
a        = 0.1    # Soft-Coulomb softening parameter
L        = 10.0   # Half-box size (a.u.)
N        = 400    # Number of grid points

max_iter  = 100   # Maximum SCF iterations
conv_tol  = 1e-6  # Convergence threshold on density-residual norm
diis_size = 6     # DIIS subspace size


# ── Helper functions ──────────────────────────────────────────────────────────

def compute_energy(eps, density, x, dx, a, R):
    """Compute the Hartree total energy.

    E = 2*eps - (1/2) * integral rho(x) V_H[rho](x) dx + V_nn

    The factor 1/2 avoids double-counting the electron-electron interaction
    that is already included in the orbital energy eps.

    Parameters
    ----------
    eps : float
        Orbital energy (lowest eigenvalue of the Hartree Fock operator).
    density : ndarray, shape (N,)
        Total electron density rho(x) = 2|phi(x)|^2.
    x : ndarray, shape (N,)
        Spatial grid.
    dx : float
        Grid spacing.
    a : float
        Soft-Coulomb softening parameter.
    R : float
        Internuclear distance (a.u.).

    Returns
    -------
    float
        Total energy in atomic units.
    """
    VH = v_hartree(density, x, dx, a)
    J  = np.trapezoid(density * VH, x)
    return 2.0 * eps - 0.5 * J + v_nn(R, a)


def run_scf(x, dx, R, a, max_iter, conv_tol, diis_size):
    """Run the Hartree SCF with DIIS density-mixing acceleration.

    The Hartree SCF equations have an unstable fixed point for the simple
    substitution iteration (due to the large self-interaction in V_H[2*phi^2]).
    DIIS (Direct Inversion of the Iterative Subspace) accelerates convergence
    by finding the optimal linear combination of past densities that minimises
    the density residual.

    Parameters
    ----------
    x : ndarray, shape (N,)
        Spatial grid.
    dx : float
        Grid spacing.
    R : float
        Internuclear distance (a.u.).
    a : float
        Soft-Coulomb softening parameter.
    max_iter : int
        Maximum number of SCF iterations.
    conv_tol : float
        Convergence threshold on the RMS density residual.
    diis_size : int
        Number of past iterates kept in the DIIS subspace.

    Returns
    -------
    phi : ndarray, shape (N,)
        Converged orbital.
    eps : float
        Converged orbital energy (a.u.).
    E : float
        Converged total energy (a.u.).
    """
    N = len(x)

    # core Hamiltonian (does not change during SCF)
    T      = kinetic_matrix(N, dx)
    H_core = T + np.diag(v_en(x, R, a))

    # initial guess: lowest eigenstate of H_core
    evals, evecs = la.eigh(H_core)
    phi = normalize(evecs[:, 0], dx)
    if phi[np.argmax(np.abs(phi))] < 0:
        phi = -phi

    density = 2.0 * phi**2     # initial density
    E_old   = 0.0

    # DIIS storage: output densities and residuals (last diis_size steps).
    # Each entry is the density PRODUCED by diagonalisation (rho_out = 2*phi²).
    # The residual is rho_in - rho_out (zero at self-consistency).
    # Extrapolation: rho_next = sum_i c_i * rho_out_i  drives residual to zero.
    densities_out = []
    residuals     = []

    print(f"{'Iter':>4}  {'eps (a.u.)':>14}  {'E_total (a.u.)':>16}  {'|res|':>12}")
    print("-" * 54)

    for it in range(max_iter):
        # ── Build Fock operator: F = H_core + V_H[rho] ───────────────────────
        # V_H uses the TOTAL density rho = 2|phi|^2, which includes
        # self-interaction (each electron sees its own field).
        VH = v_hartree(density, x, dx, a)
        F  = H_core + np.diag(VH)

        # ── Diagonalise ───────────────────────────────────────────────────────
        evals, evecs = la.eigh(F)
        eps     = evals[0]
        phi     = normalize(evecs[:, 0], dx)
        if phi[np.argmax(np.abs(phi))] < 0:
            phi = -phi

        # ── Density residual: r = rho_in - rho_out ───────────────────────────
        density_out = 2.0 * phi**2
        res         = density - density_out
        res_norm    = np.sqrt(np.trapezoid(res**2, x))

        # ── Total energy (using input density for consistency with eps) ───────
        E = compute_energy(eps, density, x, dx, a, R)
        print(f"{it + 1:>4}  {eps:>14.8f}  {E:>16.8f}  {res_norm:>12.4e}")

        if res_norm < conv_tol:
            print("-" * 54)
            print(f"SCF converged after {it + 1} iterations.")
            break

        E_old = E

        # ── DIIS update ───────────────────────────────────────────────────────
        # Store the OUTPUT density and its residual.
        densities_out.append(density_out.copy())
        residuals.append(res.copy())
        if len(densities_out) > diis_size:
            densities_out.pop(0)
            residuals.pop(0)
        m = len(densities_out)

        if m < 2:
            # Not enough history yet — use the output density directly
            # (pure substitution bootstrap for the first step).
            density = density_out
            continue

        # Build the (m+1) x (m+1) Pulay matrix B[i,j] = <r_i, r_j>
        B = np.zeros((m + 1, m + 1))
        for i in range(m):
            for j in range(m):
                B[i, j] = np.trapezoid(residuals[i] * residuals[j], x)
        B[m, :] = 1.0
        B[:, m] = 1.0
        B[m, m] = 0.0

        rhs = np.zeros(m + 1)
        rhs[m] = 1.0

        try:
            c = np.linalg.solve(B, rhs)
        except np.linalg.LinAlgError:
            # Fallback: simple average if DIIS matrix is singular
            c = np.ones(m + 1) / m
            c[m] = 0.0

        # Optimal density for next iteration: weighted combination of outputs
        density = sum(c[i] * densities_out[i] for i in range(m))

    else:
        print("-" * 54)
        print(f"WARNING: SCF did not converge within {max_iter} iterations.")

    return phi, eps, E


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    x, dx = make_grid(L, N)

    print("=" * 54)
    print("Hartree (Mean-Field) SCF for 1D H2")
    print("=" * 54)

    phi, eps, E = run_scf(x, dx, R, a, max_iter, conv_tol, diis_size)

    Vnn = v_nn(R, a)
    density = 2.0 * phi**2
    VH = v_hartree(density, x, dx, a)
    J11 = np.trapezoid(density * VH, x) / 4.0   # J_11 = (1/4) * integral rho * V_H[rho] dx

    print()
    print(f"  R                   = {R:.4f} a.u.")
    print(f"  Orbital energy  eps = {eps:.6f} a.u.")
    print(f"  V_nn                = {Vnn:.6f} a.u.")
    print(f"  J_11 (Coulomb)      = {J11:.6f} a.u.")
    print(f"  Total energy    E   = {E:.6f} a.u.")
    print(f"  (E includes self-interaction error of ~J_11 vs HF)")
    print("=" * 54)

    title = (f"Hartree Mean-Field SCF  |  R = {R} a.u.\n"
             f"E = 2*eps - J + V_nn = {E:.4f} a.u.")
    plot_orbital_and_density(x, phi, density, R, title,
                             filename="h2_mean_field.png")


if __name__ == "__main__":
    main()
