# Week 2 — Crystal Structure Prediction with Quantum ESPRESSO

This activity introduces Density Functional Theory (DFT) calculations on bulk Germanium (Ge) in the diamond cubic structure using [Quantum ESPRESSO](https://www.quantum-espresso.org/) and the PBE exchange-correlation functional. The goal is to predict the equilibrium crystal structure of Ge from first principles and compare it against experimental measurements and a database reference from the [Materials Project](https://materialsproject.org/).

## System

**Material:** Bulk Germanium — diamond cubic structure (face-centered cubic Bravais lattice, space group Fd3m)  
**Pseudopotential:** Norm-conserving pseudopotential (`Ge.upf`, UPF format)  
**Functional:** PBE (Perdew-Burke-Ernzerhof, GGA)  
**k-point sampling:** 8×8×8 Monkhorst-Pack grid  
**Wavefunction cutoff:** 60 Ry (120 Ry for vc-relax)

## Directory Structure

```
week2/
├── 00_scf/                         # Baseline SCF calculation
├── 01_make_input_experiment/       # Input from experimental structure
├── 02_make_input_materials_project/# Input from Materials Project structure
└── 03_cell_relaxation/             # Crystal structure prediction
    └── change_alat/                # Manual lattice parameter scan
```

## Activities

### `00_scf/` — Baseline SCF Calculation

A self-consistent field (SCF) calculation for bulk Ge using the experimental lattice parameter (`a = 5.658 Å`). The structure is specified compactly via `ibrav = 2` (FCC) and fractional coordinates. This serves as a reference point and introduces the basic structure of a Quantum ESPRESSO `pw.x` input file.

Key files:
- `scf.in` — PWscf input file
- `Ge.upf` — Ge pseudopotential
- `sub.cmd` — SLURM submission script (4 MPI tasks, 30 min, `instruction` partition)

### `01_make_input_experiment/` — Input from Experimental Structure

An SCF input file built from the experimentally measured crystal structure of Ge. The cell is specified explicitly via `CELL_PARAMETERS angstrom` (with `ibrav = 0`), reflecting how real structural data (e.g., from X-ray diffraction) would be used to set up a calculation. The lattice vectors correspond to the FCC primitive cell with `a ≈ 5.658 Å`.

### `02_make_input_materials_project/` — Input from Materials Project Structure

An SCF input file built from the DFT-optimized structure retrieved from the Materials Project database. As with `01_make_input_experiment`, the cell uses explicit `CELL_PARAMETERS angstrom` and `ibrav = 0`. Comparing this structure with the experimental one illustrates the typical ~1–2% overestimation of lattice parameters by GGA-PBE.

### `03_cell_relaxation/` — Crystal Structure Prediction via Variable-Cell Relaxation

A variable-cell relaxation (`vc-relax`) calculation that simultaneously optimizes atomic positions and lattice vectors using the BFGS algorithm. Starting from an initial guess of `a = 5.8 Å`, Quantum ESPRESSO minimizes the total energy with respect to all degrees of freedom at zero pressure, predicting the equilibrium crystal structure from first principles.

Key settings:
- `&IONS`: `ion_dynamics = 'bfgs'`
- `&CELL`: `cell_dynamics = 'bfgs'`, `press = 0.0 kbar`, `cell_dofree = 'all'`
- Higher wavefunction cutoff (`ecutwfc = 120 Ry`) for converged stress tensor

#### `change_alat/` — Manual Equation of State Scan

An alternative approach to finding the equilibrium lattice parameter by scanning total energies across a range of lattice constants (5.2–6.6 Å). A template input (`scf_t.in`) contains a `CHANGEME` placeholder for the lattice parameter `a`. The submission script (`sub.cmd`) loops over lattice parameter values, substituting them with `sed`, running a separate SCF calculation for each, and saving the output. The `plot.sh` script then extracts the total energies and plots Energy vs. `a` to identify the minimum, which corresponds to the predicted equilibrium lattice parameter.

## Running the Calculations

Jobs are submitted to a SLURM cluster via:

```bash
sbatch sub.cmd
```

Make sure to:
1. Replace `<cruzid>` in `sub.cmd` with your actual UCSC CruzID.
2. Create output directories (e.g., `./Ge_out/`) before running.
3. Ensure the `Ge.upf` pseudopotential is present in the working directory (or update `pseudo_dir`).

## Learning Objectives

- Write and interpret Quantum ESPRESSO `pw.x` input files for periodic solids.
- Understand the difference between `ibrav`-based and explicit (`ibrav = 0`) structural specifications.
- Compare DFT-predicted structures with experimental data and database references.
- Predict equilibrium crystal structures using variable-cell relaxation (`vc-relax`).
- Construct an equation of state by scanning total energy vs. lattice parameter.
