# Week 7 — Vibrational frequencies and phonons

## Overview

This folder contains two complete worked examples of how to compute vibrational
frequencies from first principles using **Quantum ESPRESSO** and the
**finite-difference** approach:

1. **`00_water_dimer/`** — Normal modes of a (water)$_2$ cluster in the gas phase.
   Atoms are displaced one at a time, forces are read from each SCF run, and the
   dynamical matrix is assembled and diagonalized "by hand" with custom Python
   scripts. This is the pedagogical version that mirrors the derivation in the
   lecture notes.
2. **`01_ge_phonons/`** — Phonon dispersion of bulk Ge (diamond structure) using
   **phonopy** to drive the supercell + finite-displacement workflow, and QE
   (`pw.x`) as the force calculator.

Both examples follow the same physical recipe:

1. Relax the structure so that all forces are below a tolerance.
2. Displace each atom by a small amount $\pm d$ along each Cartesian direction.
3. Run an SCF calculation for every displaced structure and collect the forces.
4. Build the dynamical matrix from a central finite difference of the forces:

   $$
   D_{Ii,Jj} = -\frac{1}{\sqrt{M_I M_J}}\,
   \frac{F_{Jj}(R_{Ii}+d/2) - F_{Jj}(R_{Ii}-d/2)}{d}
   $$

5. Diagonalize $\mathbb{D}$. Its eigenvalues are $\omega^2$ and the eigenvectors
   are the (mass-weighted) normal modes.

---

## Folder structure

```
week7/
├── 00_water_dimer/
│   ├── 00_relax/                       # Step 1: geometry relaxation of the dimer
│   │   ├── water_dimer_relax.in        # QE input (relax, 12 Å cubic box, Γ-point)
│   │   ├── O.upf, H.upf                # pseudopotentials
│   │   ├── sub.cmd                      # SLURM submission script
│   │   └── results/water_dimer_relax.out
│   └── 01_finite_displacement/         # Step 2: build & diagonalize D
│       ├── finite_displacement.py      # Generates the 36 displaced QE inputs
│       ├── compute_dynamical_matrix.py # Reads forces → D matrix → eigenvalues
│       ├── README.md                    # Quick-reference for this sub-folder
│       └── results/                     # 36 pre-computed .in/.out pairs
└── 01_ge_phonons/                      # Phonon spectrum of Ge (diamond)
    ├── ge_diamond.in                    # 8-atom conventional cell of Ge
    ├── header.in                        # &CONTROL/&SYSTEM block for supercell jobs
    ├── band.conf                        # Phonopy band-structure configuration
    ├── Ge.upf                           # Pseudopotential
    ├── README.md                        # 2×2×2 supercell recipe
    ├── results/                         # Pre-computed supercell outputs (2×2×2)
    └── 3x3x3/                           # 3×3×3 supercell variant
        ├── ge_diamond.in, header.in, band.conf
        ├── supercell.in, supercell-001.in
        ├── README.md
        └── results/                     # Pre-computed FORCE_SETS, band.yaml, etc.
```

---

## Example 1 — Water dimer (`00_water_dimer/`)

This example walks through the full algorithm explicitly: each finite displacement
is set up as an independent QE input file, and the dynamical matrix is built and
diagonalized by Python scripts you can read line by line.

### Step 1 — Relax the geometry

```bash
cd 00_water_dimer/00_relax
sbatch sub.cmd
```

The relax run uses `etot_conv_thr = 1e-5 Ry` and `forc_conv_thr = 1e-4 Ry/bohr`.
After it finishes, the relaxed coordinates appear at the end of
`water_dimer_relax.out`. You can verify convergence with:

```bash
grep '!    total energy' water_dimer_relax.out
grep -A 8 "Forces acting on atoms" water_dimer_relax.out | tail
```

> Pre-computed output: `results/water_dimer_relax.out`.

### Step 2 — Generate the displaced input files

```bash
cd ../01_finite_displacement
python finite_displacement.py ../00_relax/water_dimer_relax.out
```

`finite_displacement.py` reads the **last** structure from the relax output (the
relaxed geometry), and for each of the 6 atoms, each of the 3 Cartesian
directions, and each of the 2 signs ($\pm d$), it writes a fresh QE input.
With $d = 0.01$ Å this produces **36 input files** named
`atom_<I>_dir_<i>_<sign>.in`.

> **Python environment**: you need `ase` and `numpy` in your active environment
> before running the script.

### Step 3 — Run all 36 SCF jobs

Copy a submission script into this folder (any of the `sub.cmd` files in the
repository can be adapted — see *Cluster submission* below) and launch the jobs.
Each is a short single-point SCF, so they can be run sequentially in one batch
job or submitted as a small array.

A simple bash loop that runs them sequentially inside one SLURM job:

```bash
for f in atom_*_dir_*_*.in; do
    mpirun -np $SLURM_NTASKS pw.x -in $f > ${f%.in}.out
done
```

> Pre-computed `.out` files for all 36 displacements are in
> `01_finite_displacement/results/`.

### Step 4 — Build and diagonalize the dynamical matrix

```bash
python compute_dynamical_matrix.py
```

This script reads the forces from every `atom_<I>_dir_<i>_<sign>.out` file in
the current directory, assembles the $18\times 18$ dynamical matrix using the
central-difference formula, and prints both the matrix and its eigenvalues
(converted to cm$^{-1}$).

Expected output (PBE, see lecture notes for full table):

```
Eigenvalues of the Dynamical Matrix (cm^-2), Vibrational frequencies (cm^-1):
Eigenvalue, frequency 1:  12401330.541361,  3521.55
Eigenvalue, frequency 2:  13634948.950719,  3692.55
...
```

For 6 atoms there are 18 total degrees of freedom: **12 vibrational modes** plus
6 modes corresponding to rigid translations and rotations. The latter should
come out close to zero; small or slightly negative numerical values are
expected from finite-difference noise.

---

## Example 2 — Phonon dispersion of Ge (`01_ge_phonons/`)

For a periodic crystal we want the phonon **dispersion** $\omega(\bm{q})$, not
just frequencies at $\bm{q}=0$. We use **phonopy** to (i) build a supercell,
(ii) generate the symmetry-inequivalent finite displacements, (iii) read the
resulting forces, and (iv) Fourier-interpolate the force constants onto a
$\bm{q}$-path through the Brillouin zone.

### Prerequisites

```bash
pip install phonopy
```

Two supercell sizes are provided as separate recipes:

- `01_ge_phonons/`         → 2×2×2 supercell of the 8-atom conventional cell (64 atoms)
- `01_ge_phonons/3x3x3/`   → 3×3×3 supercell (216 atoms; better dispersion)

Both follow the **same workflow**; the only differences are the supercell
dimension in the phonopy call and the `nat` value in `header.in`.

### Step 1 — Generate the displaced supercell input(s)

```bash
cd 01_ge_phonons          # or cd 01_ge_phonons/3x3x3
phonopy --qe -d --dim="2 2 2" --pa=AUTO -c ge_diamond.in
```

For Ge in the diamond structure (high symmetry), this produces a **single**
`supercell-001.in` containing only the `CELL_PARAMETERS` and `ATOMIC_POSITIONS`
of the displaced supercell.

### Step 2 — Prepend the QE header

`phonopy --qe -d` does not write the `&CONTROL`/`&SYSTEM`/`&ELECTRONS` blocks.
We have provided a `header.in` template; concatenate it onto every
`supercell*.in`:

```bash
for file in supercell*.in; do
    cat header.in $file > tmp_ && mv tmp_ $file
done
```

> Make sure `nat` in `header.in` matches the number of atoms in the supercell
> (64 for 2×2×2, 216 for 3×3×3).

### Step 3 — Run the supercell SCF jobs

```bash
sbatch sub.cmd
```

Each supercell run is more expensive than the water-dimer SCFs — request enough
nodes/walltime accordingly.

> Pre-computed `supercell-001.out` is provided in each `results/` folder.

### Step 4 — Build the force constants

```bash
phonopy -f supercell-001.out
```

This produces `FORCE_SETS`, which contains the central-difference force constants
that phonopy will Fourier-transform onto a $\bm{q}$-grid.

### Step 5 — Plot the phonon dispersion

```bash
phonopy --qe -c ge_diamond.in -p band.conf
```

`band.conf` defines the $\bm{q}$-path (Γ → X → M → Γ → R), the dimension of the
supercell, and the primitive-axis transformation. The resulting plot has 6
branches (8 atoms / cell × 3 directions ÷ 4 since the primitive cell of diamond
has 2 atoms → 6 phonon branches: 3 acoustic + 3 optical).

---

## Cluster submission (SLURM)

A working `sub.cmd` is provided in `00_water_dimer/00_relax/`. For the Hummingbird
cluster at UCSC it looks like:

```bash
#!/bin/bash
#SBATCH -p instruction
#SBATCH -J <job_name>
#SBATCH --mail-user=<cruzid>@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o job%.j.out
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:30:00
#SBATCH --mem=8G

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

mpirun -np $SLURM_NTASKS pw.x < <input>.in > <input>.out
```

Adjust the walltime and memory upward when running the Ge supercell jobs
(`-t 02:00:00` and `--mem=16G` are reasonable starting points for the 3×3×3
case).

---

## Tips and common pitfalls

- **Always relax first.** The finite-difference formula assumes you are at the
  bottom of the potential energy surface. If forces on the relaxed structure are
  not small (≲ $10^{-4}$ Ry/bohr), the spurious linear term contaminates the
  curvature you are trying to extract.
- **Tight SCF convergence is essential.** Use `conv_thr = 1e-8` or tighter, since
  the force tolerance has to be much smaller than the change of force induced by
  the displacement $d$. Loose convergence corrupts the central difference.
- **Displacement size $d$.** $d = 0.01$ Å is the value used by both scripts; it
  is a good compromise between linear-response accuracy (small $d$) and
  numerical noise from incomplete SCF convergence (large $d$).
- **Negative eigenvalues** of the dynamical matrix are not always wrong: for
  finite molecules the 6 translation/rotation modes show up as nearly-zero
  eigenvalues, occasionally slightly negative. For crystals, however, large
  negative eigenvalues at finite $\bm{q}$ usually indicate a real instability of
  the structure.
- **Symmetry must be on** for phonopy (it uses it to reduce the number of
  displacements). Conversely, in `finite_displacement.py` we run every atom
  independently because we are not using symmetry reduction.
