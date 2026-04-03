"""
plot.py — Plotting utilities for the 1D H2 Schrödinger equation solvers.

Provides:
  plot_orbital_and_density   — 1×2 panel: orbital + electron density (1D methods)
  plot_exact_ground_state    — 2×3 panel: 2D probability density, energy levels,
                               marginal densities, diagonal wavefunction (exact method)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def plot_orbital_and_density(x, phi, density, R, title, filename=None):
    """Plot the lowest orbital and electron density for a 1D method.

    Parameters
    ----------
    x : ndarray, shape (N,)
        Spatial grid.
    phi : ndarray, shape (N,)
        Normalised spatial orbital phi(x).
    density : ndarray, shape (N,)
        Electron density rho(x) = 2*|phi(x)|^2.
    R : float
        Internuclear distance (a.u.), used to mark proton positions at ±R/2.
    title : str
        Overall figure title (method name and energy).
    filename : str or None
        If given, save figure to this path.  Otherwise display interactively.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fig.suptitle(title, fontsize=12, fontweight='bold')

    # --- orbital ---
    ax = axes[0]
    ax.plot(x, phi, color='steelblue', lw=2, label=r'$\phi(x)$')
    ax.axvline( R / 2.0, color='tomato', ls='--', lw=1.4, label='proton positions')
    ax.axvline(-R / 2.0, color='tomato', ls='--', lw=1.4)
    ax.axhline(0.0, color='k', lw=0.6, ls=':')
    ax.set_xlabel('x (a.u.)', fontsize=11)
    ax.set_ylabel(r'$\phi(x)$', fontsize=11)
    ax.set_title('Lowest Orbital', fontsize=11)
    ax.set_xlim(x[0], x[-1])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

    # --- electron density ---
    ax = axes[1]
    ax.plot(x, density, color='darkorchid', lw=2, label=r'$\rho(x) = 2|\phi|^2$')
    ax.fill_between(x, density, alpha=0.20, color='darkorchid')
    ax.axvline( R / 2.0, color='tomato', ls='--', lw=1.4, label='proton positions')
    ax.axvline(-R / 2.0, color='tomato', ls='--', lw=1.4)
    ax.set_xlabel('x (a.u.)', fontsize=11)
    ax.set_ylabel(r'$\rho(x)$', fontsize=11)
    ax.set_title('Electron Density', fontsize=11)
    ax.set_xlim(x[0], x[-1])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.4)

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved -> {filename}")
    else:
        plt.show()
    plt.close(fig)


def plot_exact_ground_state(x, psi2d, eigenvalues, R, filename=None):
    """2×3 summary plot for the exact 2-electron ground state.

    Panels:
      (0,0:2) 2D probability density |Psi(x1,x2)|^2
      (0,2)   Energy-level diagram
      (1,0)   Marginal density rho(x1)
      (1,1)   Marginal density rho(x2)
      (1,2)   Wavefunction along the diagonal x1=x2

    Parameters
    ----------
    x : ndarray, shape (N,)
        1D spatial grid (same for both electrons).
    psi2d : ndarray, shape (N, N)
        Normalised ground-state wavefunction Psi(x1, x2).
    eigenvalues : ndarray
        Array of computed eigenvalues (a.u.).
    R : float
        Internuclear distance (a.u.).
    filename : str or None
        If given, save to file; otherwise display interactively.
    """
    L = x[-1]
    prob2d = psi2d**2

    # marginal densities
    dx = x[1] - x[0]
    rho1 = np.trapezoid(prob2d, x, axis=1)   # integrate over x2 -> rho(x1)
    rho2 = np.trapezoid(prob2d, x, axis=0)   # integrate over x1 -> rho(x2)

    # diagonal slice
    N = len(x)
    diag_vals = np.array([psi2d[i, i] for i in range(N)])

    E0 = eigenvalues[0]
    Ha_to_eV = 27.2114

    fig = plt.figure(figsize=(14, 10))
    fig.suptitle(
        f"1D H$_2$ Exact Ground State  |  Soft Coulomb (R = {R} a.u.)\n"
        f"$E_0$ = {E0:+.4f} a.u. = {E0 * Ha_to_eV:.3f} eV",
        fontsize=13, fontweight='bold')

    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.44, wspace=0.38)

    # --- 2D probability density ---
    ax1 = fig.add_subplot(gs[0, :2])
    im = ax1.contourf(x, x, prob2d.T, levels=60, cmap='inferno')
    fig.colorbar(im, ax=ax1, label=r'$|\Psi(x_1, x_2)|^2$')
    ax1.axvline( R / 2.0, color='cyan', ls='--', lw=1.2, label='proton positions')
    ax1.axvline(-R / 2.0, color='cyan', ls='--', lw=1.2)
    ax1.axhline( R / 2.0, color='cyan', ls='--', lw=1.2)
    ax1.axhline(-R / 2.0, color='cyan', ls='--', lw=1.2)
    ax1.set_xlabel(r'$x_1$ (a.u.)', fontsize=11)
    ax1.set_ylabel(r'$x_2$ (a.u.)', fontsize=11)
    ax1.set_title(r'Ground-State Probability Density $|\Psi(x_1,x_2)|^2$', fontsize=11)
    ax1.legend(fontsize=9)
    ax1.set_xlim(-L, L)
    ax1.set_ylim(-L, L)

    # --- energy-level diagram ---
    ax2 = fig.add_subplot(gs[0, 2])
    for i, E in enumerate(eigenvalues):
        lw  = 2.5 if i == 0 else 1.2
        col = '#FF4500' if i == 0 else '#4477AA'
        ax2.hlines(E, 0.2, 0.8, colors=col, linewidths=lw)
        ax2.text(0.85, E, f'$E_{i}$={E:.3f}', va='center', fontsize=8, color=col)
    ax2.set_xlim(0, 1.5)
    ax2.set_xticks([])
    ax2.set_ylabel('Energy (a.u.)', fontsize=11)
    ax2.set_title('Energy Levels', fontsize=11)
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.tick_right()

    # --- marginal density rho(x1) ---
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(x, rho1, color='steelblue', lw=2)
    ax3.fill_between(x, rho1, alpha=0.25, color='steelblue')
    ax3.axvline( R / 2.0, color='tomato', ls='--', lw=1.4, label='protons')
    ax3.axvline(-R / 2.0, color='tomato', ls='--', lw=1.4)
    ax3.set_xlabel(r'$x_1$ (a.u.)', fontsize=11)
    ax3.set_ylabel(r'$\rho(x_1)$', fontsize=11)
    ax3.set_title(r'Marginal Density $\rho(x_1)$', fontsize=11)
    ax3.legend(fontsize=9)
    ax3.set_xlim(-L, L)

    # --- marginal density rho(x2) ---
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(x, rho2, color='darkorchid', lw=2)
    ax4.fill_between(x, rho2, alpha=0.25, color='darkorchid')
    ax4.axvline( R / 2.0, color='tomato', ls='--', lw=1.4, label='protons')
    ax4.axvline(-R / 2.0, color='tomato', ls='--', lw=1.4)
    ax4.set_xlabel(r'$x_2$ (a.u.)', fontsize=11)
    ax4.set_ylabel(r'$\rho(x_2)$', fontsize=11)
    ax4.set_title(r'Marginal Density $\rho(x_2)$', fontsize=11)
    ax4.legend(fontsize=9)
    ax4.set_xlim(-L, L)

    # --- diagonal wavefunction ---
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.plot(x, diag_vals, color='seagreen', lw=2)
    ax5.axhline(0, color='k', lw=0.7, ls=':')
    ax5.set_xlabel(r'$x_1 = x_2$ (a.u.)', fontsize=11)
    ax5.set_ylabel(r'$\Psi(x, x)$', fontsize=11)
    ax5.set_title(r'Wavefunction on Diagonal $x_1{=}x_2$', fontsize=11)
    ax5.set_xlim(-L, L)

    if filename is not None:
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved -> {filename}")
    else:
        plt.show()
    plt.close(fig)
