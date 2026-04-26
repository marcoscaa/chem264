# Week 5 – Equilibrium Structure of Bulk Materials and Their Surfaces

This folder contains Quantum ESPRESSO (QE) input/output files and analysis scripts for three hands-on activities covering cohesive energies, structural phase stability, and surface relaxation of crystalline Germanium (Ge).

All calculations use the PBE exchange-correlation functional with a norm-conserving pseudopotential (`Ge.upf`) and a plane-wave basis set with an 80 Ry wavefunction cutoff (60 Ry for surface calculations).

---

## Activity overview

| Folder | Topic |
|--------|-------|
| `00_cohesive_energy` | Cohesive energy of diamond and beta-Sn Ge |
| `01_change_alat` | E(V) curve of diamond Ge |
| `02_change_alat_beta_tin` | E(V) curve of beta-Sn Ge and phase diagram at T=0 |
| `03_relax_ge_diamond_surfaces` | Surface relaxation of low-index Ge(diamond) surfaces |

---

## 00_cohesive_energy — Cohesive energy of Ge

**Goal:** compute the cohesive energy of Ge in the diamond and beta-Sn phases.

The cohesive energy per atom is defined as:

$$E_\text{coh} = \frac{E_\text{bulk}}{N} - E_\text{atom}$$

where $E_\text{bulk}$ is the DFT total energy of the relaxed unit cell, $N$ is the number of atoms per cell, and $E_\text{atom}$ is the energy of an isolated atom in vacuum.

### Input files

| File | Description |
|------|-------------|
| `vc-relax_bulk_diamond.in` | Variable-cell relaxation of diamond Ge (FCC Bravais lattice, `ibrav=2`, 2 atoms/cell, 8×8×8 k-mesh) |
| `vc-relax_bulk_b-tin.in` | Variable-cell relaxation of beta-Sn Ge (BCT Bravais lattice, `ibrav=7`, 2 atoms/cell, 6×6×6 k-mesh) |
| `ge_isolated.in` | SCF calculation of a single Ge atom in a 10 Å cubic box. Spin polarization (`nspin=2`) is enabled to correctly handle the unpaired electrons of the isolated atom. Only the Γ point is needed. |

### Key concepts

- The `vc-relax` calculation simultaneously relaxes ionic positions and the unit-cell lattice parameters until the residual forces and stress are below the convergence thresholds.
- The isolated-atom calculation uses a large supercell to avoid spurious self-interaction between periodic images.
- Spin polarization must be activated for open-shell systems (isolated atoms with unpaired electrons).
- The experimental cohesive energy of Ge (diamond) is −3.75 eV/atom; PBE is expected to slightly underestimate this.

---

## 01_change_alat — Energy–volume curve of diamond Ge

**Goal:** map the total energy as a function of the lattice constant `a` for the diamond structure of Ge, and identify the equilibrium lattice constant as the minimum of E(a).

### Files

| File | Description |
|------|-------------|
| `scf_t.in` | Template SCF input file. The placeholder `CHANGEME` is replaced by the actual value of `a`. Atomic coordinates are given in `alat` units so that the entire structure scales uniformly with `a`. |
| `sub.cmd` | SLURM script that loops over 21 values of `a` (5.2–6.6 Å), generates a complete input file via `sed`, runs `pw.x`, and saves one output file per `a`. |
| `plot.sh` | Shell script that extracts the converged total energy from each output file and generates `energy_alat.pdf` using matplotlib. |
| `results/` | Pre-computed output files (`scf_<alat>.out`) for each lattice constant. |

### Key concepts

- Writing atomic positions in `alat` units is essential: as `a` changes, the fractional coordinates remain fixed and the Cartesian positions scale proportionally, keeping the crystal geometry intact.
- The energy minimum determines the DFT equilibrium lattice constant without a full `vc-relax`.
- Energy extraction from QE output files uses the line starting with `!` (the converged SCF energy).

---

## 02_change_alat_beta_tin — Phase stability and T=0 phase diagram

**Goal:** compute E(V) curves for both the diamond and beta-Sn phases of Ge, then determine the thermodynamically stable phase as a function of applied hydrostatic pressure at T=0 using the enthalpy:

$$H(P) = E + PV$$

### Files

| File | Description |
|------|-------------|
| `scf_t.in` | Template SCF input for beta-Sn Ge (BCT cell, `ibrav=7`). The placeholder `CHANGEME` is the `celldm(1)` parameter (lattice constant in Bohr). The c/a ratio is held fixed at 0.537. |
| `vc-relax_bulk.in` | Variable-cell relaxation of beta-Sn Ge to find the equilibrium geometry. |
| `sub.cmd` | SLURM script that loops over values of `a`, converts Å to Bohr using Python, and runs QE for each volume. |
| `plot.sh` | Generates two plots: (1) cohesive energy E(a) for both phases; (2) enthalpy H(a,P) for both phases at a user-specified pressure P (in GPa). The crossing of the two enthalpy curves gives the predicted phase-transition pressure. |
| `results/` | Pre-computed output files. |

### Key concepts

- The diamond and beta-Sn phases have different densities (the beta-Sn phase is denser). Under high pressure, the PV term disfavors the lower-density diamond phase.
- The enthalpy is minimized at fixed P to find the equilibrium volume for each phase; the phase with the lower enthalpy minimum is the stable one.
- The phase-transition pressure is where $H_\text{diamond}(P) = H_\text{beta-Sn}(P)$. DFT/PBE predicts a transition near ~10 GPa, consistent with experiment.
- The top-level `plot.sh` script (in this folder) collects data from both `01_change_alat` and `02_change_alat_beta_tin` and produces the combined energy and enthalpy plots.

---

## 03_relax_ge_diamond_surfaces — Surface relaxation of Ge

**Goal:** study the atomic-scale relaxation and reconstruction of low-index surfaces of diamond Ge using DFT `relax` calculations (fixed cell, mobile ions). Surface energies are used to assess convergence with respect to the surface supercell size.

### Sub-folders

```
03_relax_ge_diamond_surfaces/
├── 100/           # (100) surface, multiple supercell sizes
│   ├── 1x1/       # 1×1 surface cell (8 atoms, 4×4×1 k-mesh)
│   ├── 1x2/       # 1×2 supercell
│   ├── 2x1/       # 2×1 supercell
│   ├── 2x2/       # 2×2 supercell
│   ├── 4x2/       # 4×2 supercell
│   └── 4x4/       # 4×4 supercell
├── 110/           # (110) surface (96 atoms)
└── 111/           # (111) surface (128 atoms, Γ-only k-sampling)
```

Each subfolder contains `surface_relax.in` (QE `relax` input) and `results/surface_relax.out`.

### Scripts

| File | Description |
|------|-------------|
| `create_pw_input_surface_from_materials_project.py` | Fetches the crystal structure from the Materials Project (using the mp-id, e.g., `mp-32` for Ge), constructs a slab with a given Miller index and vacuum (~10 Å), replicates it into a surface supercell of size nA×nB, adds small random displacements to break symmetry, and writes a QE `relax` input file. |
| `compute_surface_energy.py` | Reads the bulk and slab QE output files (via ASE) and computes the surface energy in J/m²: $\gamma = (E_\text{slab} - N \cdot \varepsilon_\text{bulk}) / (2A)$ |

**Usage:**
```bash
# Create input file for the (100) 2×1 surface
python create_pw_input_surface_from_materials_project.py mp-32 100 2 1

# Compute surface energy
python compute_surface_energy.py ../00_cohesive_energy/results/vc-relax_bulk_diamond.out results/surface_relax.out
```

### Key concepts

- **Slab model:** a finite-thickness slab (~4 layers, ~10 Å vacuum) is used to simulate the surface while preserving 3D periodic boundary conditions. The slab must be thick enough for the central layers to be bulk-like.
- **Surface relaxation:** after geometry optimization, atoms near the surface move from their bulk positions to minimize dangling bonds. In larger supercells (e.g., 2×2 or larger), atoms can undergo **surface reconstruction**, involving significant lateral displacements and dimerization.
- **Supercell convergence:** the surface energy must converge with respect to both slab thickness and surface supercell area. The (100) surface of Ge is known to reconstruct with a 2×1 periodicity (dimer rows).
- **Surface energy formula:** the factor of ½ accounts for the two surfaces present in a periodic slab calculation.
- The k-point mesh is dense along the surface directions (x, y) but only 1 point along z (normal to the surface), reflecting the broken periodicity in that direction.
