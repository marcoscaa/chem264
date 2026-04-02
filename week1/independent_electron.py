"""
independent_electron.py — Independent-electron approximation for 1D H2.

Approximation
-------------
Each electron moves independently in the field of the two protons.
Electron-electron repulsion is completely neglected.  The two electrons
occupy the same lowest-energy spatial orbital phi_0(x) (closed-shell
singlet), giving a total energy

    E = 2 * eps_0 + V_nn

where eps_0 is the lowest eigenvalue of the single-electron Hamiltonian

    h = T + V_en(x)

and V_nn is the nuclear-nuclear repulsion (a constant added at the end).

This is the crudest level of theory: it overestimates binding because it
omits the destabilising electron-electron Coulomb repulsion.

Usage
-----
    python independent_electron.py
"""

import numpy as np
import scipy.linalg as la

from potentials import v_en, v_nn
from utils import make_grid, kinetic_matrix, normalize
from plot import plot_orbital_and_density

# ── Parameters ────────────────────────────────────────────────────────────────
R   = 1.4    # Internuclear distance (a.u.)
a   = 0.1    # Soft-Coulomb softening parameter
L   = 10.0   # Half-box size (a.u.)
N   = 400    # Number of grid points


# ── Helper functions ──────────────────────────────────────────────────────────

def build_hamiltonian(x, dx, R, a):
    """Build the single-electron Hamiltonian matrix h = T + diag(V_en).

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

    Returns
    -------
    H : ndarray, shape (N, N)
        Single-electron Hamiltonian matrix.
    """
    T = kinetic_matrix(N, dx)
    V = v_en(x, R, a)
    H = T + np.diag(V)
    return H


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # --- grid ---
    x, dx = make_grid(L, N)

    # --- Hamiltonian ---
    H = build_hamiltonian(x, dx, R, a)

    # --- diagonalise (eigh exploits symmetry, returns ascending eigenvalues) ---
    eigenvalues, eigenvectors = la.eigh(H)

    # --- ground-state orbital ---
    phi = eigenvectors[:, 0]
    phi = normalize(phi, dx)

    # enforce positive phase convention
    if phi[np.argmax(np.abs(phi))] < 0:
        phi = -phi

    eps0 = eigenvalues[0]

    # --- total energy ---
    Vnn   = v_nn(R, a)
    E_tot = 2.0 * eps0 + Vnn

    # --- print results ---
    print("=" * 50)
    print("Independent-Electron Approximation for 1D H2")
    print("=" * 50)
    print(f"  R            = {R:.4f} a.u.")
    print(f"  Orbital energy eps_0 = {eps0:.6f} a.u.")
    print(f"  V_nn         = {Vnn:.6f} a.u.")
    print(f"  Total energy E = 2*eps_0 + V_nn = {E_tot:.6f} a.u.")
    print("=" * 50)

    # --- plot ---
    density = 2.0 * phi**2
    title = (f"Independent Electron  |  R = {R} a.u.\n"
             f"E = 2*eps_0 + V_nn = {E_tot:.4f} a.u.")
    plot_orbital_and_density(x, phi, density, R, title,
                             filename="h2_independent_electron.png")


if __name__ == "__main__":
    main()
