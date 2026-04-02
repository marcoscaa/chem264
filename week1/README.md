# 1D H2 Molecule — Schrodinger Equation Solvers

A set of Python scripts that solve the electronic Schrodinger equation for a
one-dimensional model of the H2 molecule at four levels of theory, from the
crudest independent-electron picture to the numerically exact two-electron
solution.

## Model

All scripts share the same physical model (atomic units: hbar = m_e = e = 1):

- **Born-Oppenheimer approximation** — the two protons are fixed at positions
  x = ±R/2 with R = 1.4 a.u. (near the H2 equilibrium bond length).
- **Soft Coulomb potential** — the singular 1/r interaction is regularised as
  `V(r) = 1/sqrt(r^2 + a^2)` with softening parameter a = 0.5 a.u.
  This avoids the 1D collapse problem while preserving the qualitative physics.
- **Box** — electrons are confined to a box of half-width L = 10 a.u.
- **Finite-difference grid** — 1D scripts use N = 400 points; the exact 2D
  solver uses N = 128 points per dimension.

## Methods

| Script | Method | Description |
|---|---|---|
| `independent_electron.py` | Independent electron | No electron-electron interaction. Each electron sees only the nuclear potential. Total energy E = 2*eps_0 + V_nn. |
| `mean_field.py` | Hartree SCF | Electrons interact through the mean Hartree potential V_H[rho]. SCF uses orbital (density) mixing. |
| `hartree_fock.py` | Restricted Hartree-Fock | Formally identical to Hartree for a single doubly-occupied real orbital (K = J), but uses Fock-matrix mixing for better convergence. Reports J and K explicitly. |
| `numerical_exact.py` | Full 2-electron exact | Both electrons treated simultaneously on a 2D grid. Includes all electron correlation. Sparse ARPACK diagonalisation. |

Shared utilities:

| Script | Contents |
|---|---|
| `potentials.py` | `soft_coulomb`, `v_en`, `v_nn`, `v_hartree` |
| `utils.py` | `make_grid`, `kinetic_matrix`, `normalize` |
| `plot.py` | `plot_orbital_and_density`, `plot_exact_ground_state` |

## Requirements

```
numpy
scipy
matplotlib
```

Install with:

```bash
pip install numpy scipy matplotlib
```

## How to Run

Each script is self-contained and runs from the command line:

```bash
python independent_electron.py
python mean_field.py
python hartree_fock.py
python numerical_exact.py   # slowest (~30 s for N=128)
```

## Output Files

| File | Produced by | Contents |
|---|---|---|
| `h2_independent_electron.png` | `independent_electron.py` | Orbital and electron density |
| `h2_mean_field.png` | `mean_field.py` | Converged Hartree orbital and density |
| `h2_hartree_fock.png` | `hartree_fock.py` | Converged RHF orbital and density |
| `h2_exact.png` | `numerical_exact.py` | 2D probability density, energy levels, marginal densities, diagonal wavefunction |
| `eigenvalues.txt` | `numerical_exact.py` | Table of the lowest 6 two-electron eigenvalues in a.u. and eV |

- **Independent electron** (most negative) — completely ignores electron-electron
  repulsion, so it vastly overestimates binding.
- **Exact** — treats both electrons simultaneously on a 2D grid, recovering all
  electron correlation. This is the true benchmark.
- **Hartree-Fock** — includes the average electron-electron repulsion (Coulomb
  integral J₁₁) and the exchange interaction (K₁₁), which exactly cancels the
  self-interaction error present in the Hartree method.
- **Hartree** (least negative, highest energy) — includes electron-electron
  repulsion via the mean Hartree potential V_H[ρ] built from the *total* density
  ρ = 2|φ|². This means each electron interacts with itself (self-interaction
  error), inflating the energy by an extra J₁₁ compared to HF:
  Energy: E_H = 2·h₁₁ + 2·J₁₁ + V_nn = E_HF + J₁₁ > E_HF

