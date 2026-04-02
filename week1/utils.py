"""
utils.py — Grid and linear-algebra utilities for the 1D H2 solvers.

Provides:
  make_grid     — uniform 1D grid on [-L, L)
  kinetic_matrix — finite-difference kinetic energy matrix (-1/2 d^2/dx^2)
  normalize      — L2-normalise a wavefunction on the grid
"""

import numpy as np


def make_grid(L, N):
    """Create a uniform spatial grid on [-L, L).

    Uses endpoint=False so that periodic copies of the box edge do not
    overlap (consistent with a particle-in-a-box boundary condition where
    the wavefunction vanishes at ±L).

    Parameters
    ----------
    L : float
        Half-box size (a.u.).
    N : int
        Number of grid points.

    Returns
    -------
    x : ndarray, shape (N,)
        Grid points.
    dx : float
        Uniform grid spacing.
    """
    x = np.linspace(-L, L, N, endpoint=False)
    dx = x[1] - x[0]
    return x, dx


def kinetic_matrix(N, dx):
    """Dense finite-difference matrix for T = -1/2 d^2/dx^2.

    Uses the standard three-point central-difference formula:
        T[i, i]   =  1 / dx^2
        T[i, i±1] = -0.5 / dx^2

    Parameters
    ----------
    N : int
        Number of grid points.
    dx : float
        Grid spacing.

    Returns
    -------
    T : ndarray, shape (N, N)
        Kinetic energy matrix in atomic units.
    """
    diag     = np.ones(N) / dx**2
    off_diag = -0.5 * np.ones(N - 1) / dx**2
    T = np.diag(diag) + np.diag(off_diag, k=1) + np.diag(off_diag, k=-1)
    return T


def normalize(psi, dx):
    """Return psi normalised so that sum(psi^2)*dx = 1.

    Parameters
    ----------
    psi : ndarray, shape (N,)
        Wavefunction (real).
    dx : float
        Grid spacing.

    Returns
    -------
    ndarray, shape (N,)
        Normalised wavefunction.
    """
    norm = np.sqrt(np.sum(psi**2) * dx)
    return psi / norm
