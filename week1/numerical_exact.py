"""
numerical_exact.py — Exact 2-electron solution for 1D H2.

Model
-----
Both electrons (coordinates x1, x2) are treated explicitly on a 2D grid.
The full two-electron Hamiltonian in atomic units (hbar = m_e = e = 1) is

    H = -1/2 d^2/dx1^2  -  1/2 d^2/dx2^2
        + V_en(x1) + V_en(x2)          (electron-nucleus attraction)
        + V_ee(x1, x2)                 (electron-electron repulsion)
        + V_nn                         (nucleus-nucleus repulsion, constant)

with soft Coulomb interactions:

    V_en(xi)     = -1/sqrt((xi - R/2)^2 + a^2)  -  1/sqrt((xi + R/2)^2 + a^2)
    V_ee(x1,x2)  = +1/sqrt((x1 - x2)^2 + a^2)
    V_nn         = +1/sqrt(R^2 + a^2)

The N^2 x N^2 Hamiltonian is assembled using Kronecker products of the
sparse 1D kinetic matrix plus a diagonal potential matrix.  The lowest
n_eigs eigenstates are computed with scipy's ARPACK wrapper (eigsh).

This solution is numerically exact within the finite-difference
discretisation and the soft-Coulomb model: it includes all electron
correlation effects.

Usage
-----
    python numerical_exact.py
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from potentials import v_en, v_nn, soft_coulomb
from utils import make_grid
from plot import plot_exact_ground_state

# ── Parameters ────────────────────────────────────────────────────────────────
R      = 1.4    # Internuclear distance (a.u.)
a      = 0.1    # Soft-Coulomb softening parameter
L      = 10.0   # Half-box size (a.u.)
N      = 128    # Grid points per dimension  (Hamiltonian is N^2 x N^2)
n_eigs = 6      # Number of lowest eigenstates to compute


# ── Helper functions ──────────────────────────────────────────────────────────

def build_kinetic_1d_sparse(N, dx):
    """Sparse finite-difference matrix for T = -1/2 d^2/dx^2 in 1D.

    Parameters
    ----------
    N : int
        Number of grid points.
    dx : float
        Grid spacing.

    Returns
    -------
    scipy.sparse.csr_matrix, shape (N, N)
        Sparse kinetic energy matrix.
    """
    diag     =  np.ones(N)     / dx**2
    off_diag = -0.5 * np.ones(N - 1) / dx**2
    T1d = sp.diags([off_diag, diag, off_diag], [-1, 0, 1],
                   shape=(N, N), format='csr')
    return T1d


def build_hamiltonian(x, dx, N, R, a):
    """Assemble the full 2-electron sparse Hamiltonian.

    Parameters
    ----------
    x : ndarray, shape (N,)
        1D spatial grid.
    dx : float
        Grid spacing.
    N : int
        Number of grid points per dimension.
    R : float
        Internuclear distance (a.u.).
    a : float
        Soft-Coulomb softening parameter.

    Returns
    -------
    H : scipy.sparse.csr_matrix, shape (N^2, N^2)
        Full 2-electron Hamiltonian.
    """
    # 1D sparse kinetic matrix
    T1d = build_kinetic_1d_sparse(N, dx)

    # 2-electron kinetic energy via Kronecker products
    I   = sp.eye(N, format='csr')
    T2e = sp.kron(T1d, I) + sp.kron(I, T1d)

    # 2D potential on the N x N mesh
    X1, X2 = np.meshgrid(x, x, indexing='ij')  # X1[i,j]=x[i], X2[i,j]=x[j]

    # electron-nuclear attraction (uses potentials.v_en for each coordinate)
    Ven2d = v_en(X1, R, a) + v_en(X2, R, a)

    # electron-electron repulsion
    Vee2d = soft_coulomb(X1 - X2, a)

    # nuclear-nuclear repulsion (constant)
    Vnn   = v_nn(R, a)

    # total potential flattened to a vector
    V_vec = (Ven2d + Vee2d + Vnn).ravel()
    V_mat = sp.diags(V_vec, 0, format='csr')

    H = T2e + V_mat
    return H


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    x, dx = make_grid(L, N)

    print("=" * 60)
    print("Exact 2-Electron Solution for 1D H2")
    print("=" * 60)
    print(f"  Grid: {N} x {N} points, dx = {dx:.4f} a.u., box = [{-L}, {L}] a.u.")

    H = build_hamiltonian(x, dx, N, R, a)
    print(f"  Hamiltonian size: {H.shape[0]} x {H.shape[1]},  "
          f"non-zeros: {H.nnz:,}")

    print(f"\n  Solving for {n_eigs} lowest eigenvalues via ARPACK ...")
    eigenvalues, eigenvectors = spla.eigsh(H, k=n_eigs, which='SA')

    # sort ascending (eigsh usually does, but guarantee it)
    idx          = np.argsort(eigenvalues)
    eigenvalues  = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # --- print eigenvalues ---
    Ha_to_eV = 27.2114
    print()
    print(f"{'Index':>6}  {'Energy (a.u.)':>18}  {'Energy (eV)':>14}")
    print("-" * 44)
    for i, E in enumerate(eigenvalues):
        print(f"{i:>6}  {E:>+18.6f}  {E * Ha_to_eV:>14.4f}")
    print()
    print(f"  Ground-state energy E0 = {eigenvalues[0]:+.6f} a.u. "
          f"= {eigenvalues[0] * Ha_to_eV:.4f} eV")
    print("=" * 60)

    # --- save eigenvalues ---
    ev_file = "eigenvalues.txt"
    header  = (f"Eigenvalues for 1D H2 (soft Coulomb, a={a}, R={R} a.u.)\n"
               f"{'Index':>6}  {'Energy (a.u.)':>18}  {'Energy (eV)':>14}")
    np.savetxt(ev_file,
               np.column_stack([np.arange(n_eigs),
                                eigenvalues,
                                eigenvalues * Ha_to_eV]),
               header=header,
               fmt="%6d  %18.8f  %14.6f")
    print(f"  Eigenvalues saved -> {ev_file}")

    # --- ground-state wavefunction ---
    psi0 = eigenvectors[:, 0].reshape(N, N)
    # normalise: sum(psi^2) * dx^2 = 1
    psi0 = psi0 / np.sqrt(np.sum(psi0**2) * dx**2)

    # --- plot ---
    plot_exact_ground_state(x, psi0, eigenvalues, R,
                            filename="h2_exact.png")


if __name__ == "__main__":
    main()
