"""
energy_curves.py — Potential energy curves for 1D H2 at four levels of theory.

Sweeps over a range of internuclear distances R and computes the total energy
at each level:

  1. Independent electron  — no e-e repulsion; E = 2*eps_0 + V_nn
  2. Hartree (mean-field)  — Hartree SCF with DIIS density mixing
  3. Hartree-Fock (RHF)    — DIIS density mixing with the HF Fock operator
  4. Numerically exact     — full 2-electron sparse diagonalisation (ARPACK)

All methods share the same soft-Coulomb model and utility functions from
potentials.py and utils.py.

Warm starts
-----------
All three SCF methods and the exact solver are seeded with the converged
solution from the previous R value.  This keeps each solver on the same
solution branch across the scan and greatly reduces iteration counts.

Grid aliasing in the exact solver
-----------------------------------
With softening parameter a = 0.1, the nuclear potential V_en has a half-width
of ~a = 0.1 a.u.  When a nucleus coincides with a grid point the diagonal
matrix element is -1/a = -10 a.u.; halfway between grid points it drops to
roughly -6 a.u.  This ~4 a.u. oscillation corrupts the energy curve.

Fix: for each R the exact energy is averaged over n_shifts = 4 uniform grid
translations spaced dx/n_shifts apart.  Because the aliasing error is a
periodic function of (R/2 mod dx), averaging over one full period cancels the
dominant harmonic.  The ARPACK starting vector (v0) is recycled across both
grid shifts and R values, making each warm eigsh call typically 10-50x faster
than a cold start.

Grid sizes
----------
  N_scf   = 300   (for SCF methods)
  L_scf   = 12.0  (half-box size for SCF; enlarged to accommodate large R)
  N_exact = 128   (per dimension; same as numerical_exact.py)
  L_exact = 10.0  (half-box size for exact solver)

Output
------
  h2_energy_curves.png  — potential energy curves for all four methods

Usage
-----
    python energy_curves.py
"""

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

from potentials import v_en, v_nn, v_hartree, soft_coulomb
from utils import make_grid, kinetic_matrix, normalize


# ── Shared parameters ─────────────────────────────────────────────────────────
a        = 0.1    # soft-Coulomb softening parameter
L_scf    = 12.0   # half-box size for SCF methods (a.u.)
L_exact  = 10.0   # half-box size for exact solver (a.u.)
N_scf    = 300    # grid points for SCF methods
N_exact  = 128    # grid points per dimension for exact solver
n_shifts = 4      # grid translations for aliasing correction (exact only)

R_values = np.linspace(0.5, 7.0, 25)   # bond distances to scan (a.u.)


# ── Independent electron ──────────────────────────────────────────────────────

def energy_independent(x, dx, R, a):
    """Total energy in the independent-electron approximation.

    E = 2*eps_0 + V_nn  (no electron-electron repulsion).

    Returns
    -------
    float  Total energy (a.u.).
    """
    N    = len(x)
    T    = kinetic_matrix(N, dx)
    H    = T + np.diag(v_en(x, R, a))
    eps0 = la.eigh(H, eigvals_only=True)[0]
    return 2.0 * eps0 + v_nn(R, a)


# ── DIIS helper ───────────────────────────────────────────────────────────────

def _diis_update(densities_out, residuals, density_out, res):
    """Store new density/residual and return the DIIS-extrapolated density."""
    densities_out.append(density_out.copy())
    residuals.append(res.copy())

    m = len(densities_out)
    if m < 2:
        return density_out

    B = np.zeros((m + 1, m + 1))
    for i in range(m):
        for j in range(m):
            # inner product uses trapezoidal rule implied by the spacing
            B[i, j] = np.dot(residuals[i], residuals[j])
    B[m, :] = 1.0
    B[:, m] = 1.0
    B[m, m] = 0.0
    rhs = np.zeros(m + 1)
    rhs[m] = 1.0
    try:
        c = np.linalg.solve(B, rhs)
    except np.linalg.LinAlgError:
        c = np.ones(m + 1) / m
        c[m] = 0.0
    return sum(c[i] * densities_out[i] for i in range(m))


# ── Hartree SCF (DIIS) ────────────────────────────────────────────────────────

def energy_hartree(x, dx, R, a,
                   density_init=None,
                   max_iter=300, conv_tol=1e-6, diis_size=8):
    """Total Hartree energy via DIIS-accelerated density-mixing SCF.

    Fock operator: F_H = h + V_H[rho]   (rho = 2|phi|^2)
    Energy:        E   = 2*eps - (1/2)*integral(rho * V_H[rho]) + V_nn

    Parameters
    ----------
    density_init : ndarray or None
        Warm-start density (from previous R).  None uses H_core ground state.

    Returns
    -------
    density : ndarray   Converged density (pass back as density_init).
    E       : float     Total energy (a.u.), or NaN if not converged.
    """
    N      = len(x)
    T      = kinetic_matrix(N, dx)
    H_core = T + np.diag(v_en(x, R, a))

    if density_init is not None:
        density = density_init.copy()
    else:
        evals, evecs = la.eigh(H_core)
        phi     = normalize(evecs[:, 0], dx)
        if phi[np.argmax(np.abs(phi))] < 0:
            phi = -phi
        density = 2.0 * phi**2

    densities_out = []
    residuals     = []

    for _ in range(max_iter):
        VH  = v_hartree(density, x, dx, a)
        F   = H_core + np.diag(VH)
        evals, evecs = la.eigh(F)
        eps = evals[0]
        phi = normalize(evecs[:, 0], dx)
        if phi[np.argmax(np.abs(phi))] < 0:
            phi = -phi

        density_out = 2.0 * phi**2
        res         = density - density_out
        res_norm    = np.sqrt(dx * np.dot(res, res))

        if res_norm < conv_tol:
            VH = v_hartree(density_out, x, dx, a)
            J  = np.trapz(density_out * VH, x)
            return density_out, 2.0 * eps - 0.5 * J + v_nn(R, a)

        if len(densities_out) >= diis_size:
            densities_out.pop(0)
            residuals.pop(0)
        density = _diis_update(densities_out, residuals, density_out, res)

    return density_out, np.nan


# ── Hartree-Fock SCF (DIIS) ───────────────────────────────────────────────────

def energy_hf(x, dx, R, a,
              density_init=None,
              max_iter=300, conv_tol=1e-6, diis_size=8):
    """Total RHF energy via DIIS-accelerated density-mixing SCF.

    For a single doubly-occupied real orbital the HF Fock operator reduces to
        F_HF = h + V_H[phi^2] = h + (1/2) * V_H[rho]
    which differs from the Hartree operator only by the factor 1/2 in V_H.
    This factor exactly cancels the self-interaction (exchange = Coulomb for a
    single real orbital).

    Energy:  E = 2*eps - K_11 + V_nn,  K_11 = J_11 = (1/4) * integral(rho * V_H[rho])

    Parameters
    ----------
    density_init : ndarray or None
        Warm-start density.  None uses H_core ground state.

    Returns
    -------
    density : ndarray   Converged density.
    E       : float     Total energy (a.u.), or NaN if not converged.
    """
    N      = len(x)
    T      = kinetic_matrix(N, dx)
    H_core = T + np.diag(v_en(x, R, a))

    if density_init is not None:
        density = density_init.copy()
    else:
        evals, evecs = la.eigh(H_core)
        phi     = normalize(evecs[:, 0], dx)
        if phi[np.argmax(np.abs(phi))] < 0:
            phi = -phi
        density = 2.0 * phi**2

    densities_out = []
    residuals     = []

    for _ in range(max_iter):
        # HF Fock operator: h + V_H[rho/2] = h + (1/2)*V_H[rho]
        VH_half = v_hartree(density / 2.0, x, dx, a)
        F       = H_core + np.diag(VH_half)
        evals, evecs = la.eigh(F)
        eps = evals[0]
        phi = normalize(evecs[:, 0], dx)
        if phi[np.argmax(np.abs(phi))] < 0:
            phi = -phi

        density_out = 2.0 * phi**2
        res         = density - density_out
        res_norm    = np.sqrt(dx * np.dot(res, res))

        if res_norm < conv_tol:
            VH  = v_hartree(density_out, x, dx, a)
            J   = np.trapz(density_out * VH, x)   # = 4*J_11
            K11 = J / 4.0                          # = J_11 for single real orbital
            return density_out, 2.0 * eps - K11 + v_nn(R, a)

        if len(densities_out) >= diis_size:
            densities_out.pop(0)
            residuals.pop(0)
        density = _diis_update(densities_out, residuals, density_out, res)

    return density_out, np.nan


# ── Numerically exact (full 2-electron) ───────────────────────────────────────

def _kinetic_1d_sparse(N, dx):
    """Sparse 3-point finite-difference kinetic matrix for -1/2 d²/dx²."""
    diag     =  np.ones(N) / dx**2
    off_diag = -0.5 * np.ones(N - 1) / dx**2
    return sp.diags([off_diag, diag, off_diag], [-1, 0, 1],
                    shape=(N, N), format='csr')


def energy_exact(x, dx, N, R, a, n_shifts=4, v0=None):
    """Ground-state energy averaged over n_shifts uniform grid translations.

    The ARPACK starting vector v0 is recycled across grid shifts and R values.
    A warm v0 typically reduces the Lanczos iteration count by 10-50x.

    Parameters
    ----------
    v0 : ndarray or None
        Starting vector for ARPACK.  Pass the vector returned by the previous
        call to reuse across R values.

    Returns
    -------
    E   : float     Anti-aliased ground-state total energy (a.u.).
    v0  : ndarray   Last converged eigenvector (for warm start at next R).
    """
    T1d = _kinetic_1d_sparse(N, dx)
    I   = sp.eye(N, format='csr')
    T2e = sp.kron(T1d, I) + sp.kron(I, T1d)

    E_sum  = 0.0
    v0_cur = v0

    for k in range(n_shifts):
        xs = x + k * dx / n_shifts   # translate the entire grid
        X1, X2 = np.meshgrid(xs, xs, indexing='ij')
        V_vec  = (v_en(X1, R, a) + v_en(X2, R, a)
                  + soft_coulomb(X1 - X2, a)
                  + v_nn(R, a)).ravel()
        H = T2e + sp.diags(V_vec, format='csr')

        kwargs = dict(k=1, which='SA')
        if v0_cur is not None and v0_cur.shape[0] == H.shape[0]:
            kwargs['v0'] = v0_cur
        evals, evecs = spla.eigsh(H, **kwargs)

        E_sum  += evals[0]
        v0_cur  = evecs[:, 0]   # recycle for next shift / next R

    return E_sum / n_shifts, v0_cur


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    x_scf,   dx_scf   = make_grid(L_scf,   N_scf)
    x_exact, dx_exact = make_grid(L_exact, N_exact)

    n_R    = len(R_values)
    E_ind  = np.full(n_R, np.nan)
    E_hart = np.full(n_R, np.nan)
    E_hf   = np.full(n_R, np.nan)
    E_ex   = np.full(n_R, np.nan)

    print(f"Scanning {n_R} bond distances: R = {R_values[0]:.2f} … {R_values[-1]:.2f} a.u.")
    print(f"SCF grid  : N={N_scf},  L={L_scf}   (dx={dx_scf:.4f} a.u.)")
    print(f"Exact grid: {N_exact}x{N_exact}, L={L_exact}   "
          f"(dx={dx_exact:.4f} a.u., {n_shifts}-shift aliasing correction)\n")
    print(f"{'R (a.u.)':>10}  {'Indep':>12}  {'Hartree':>12}  "
          f"{'HF':>12}  {'Exact':>12}")
    print("-" * 66)

    # warm-start containers
    density_hart = None
    density_hf   = None
    v0_exact     = None

    for i, R in enumerate(R_values):
        E_ind[i] = energy_independent(x_scf, dx_scf, R, a)

        density_hart, E_hart[i] = energy_hartree(
            x_scf, dx_scf, R, a, density_init=density_hart)
        if np.isnan(E_hart[i]):          # warm start failed — try cold start
            density_hart, E_hart[i] = energy_hartree(
                x_scf, dx_scf, R, a, density_init=None)

        density_hf, E_hf[i] = energy_hf(
            x_scf, dx_scf, R, a, density_init=density_hf)
        if np.isnan(E_hf[i]):            # warm start failed — try cold start
            density_hf, E_hf[i] = energy_hf(
                x_scf, dx_scf, R, a, density_init=None)

        E_ex[i], v0_exact = energy_exact(
            x_exact, dx_exact, N_exact, R, a, n_shifts, v0=v0_exact)

        print(f"{R:10.3f}  {E_ind[i]:12.6f}  {E_hart[i]:12.6f}  "
              f"{E_hf[i]:12.6f}  {E_ex[i]:12.6f}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print()
    print("  Energies at R = 1.4 a.u. (reference geometry):")
    R_ref = 1.4
    for label, E in [("Independent", E_ind), ("Hartree", E_hart),
                     ("HF", E_hf), ("Exact", E_ex)]:
        idx = np.nanargmin(np.abs(R_values - R_ref))
        if not np.isnan(E[idx]):
            print(f"    {label:<14} E(R~{R_values[idx]:.2f}) = {E[idx]:.6f} a.u.")

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))

    styles = [
        ("Independent electron", E_ind,  "tab:blue",   "--",  1.5),
        ("Hartree (mean-field)", E_hart, "tab:orange",  "-.",  1.5),
        ("Hartree-Fock",         E_hf,  "tab:green",   ":",   1.8),
        ("Exact",                E_ex,  "tab:red",     "-",   2.0),
    ]
    for label, E, color, ls, lw in styles:
        mask = ~np.isnan(E)
        ax.plot(R_values[mask], E[mask],
                color=color, linestyle=ls, linewidth=lw, label=label)

    ax.set_xlabel("Internuclear distance $R$ (a.u.)", fontsize=12)
    ax.set_ylabel("Total energy (a.u.)", fontsize=12)
    ax.set_title("Potential energy curves for 1D H$_2$\n"
                 f"(soft Coulomb, $a$ = {a} a.u.)", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    fname = "h2_energy_curves.png"
    plt.savefig(fname, dpi=150)
    print(f"\nPlot saved -> {fname}")
    plt.show()


if __name__ == "__main__":
    main()
