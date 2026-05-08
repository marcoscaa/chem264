# Week 6 — Elastic Constants of Germanium via Finite Differences

## Overview

This folder contains Quantum ESPRESSO (PW) input files to compute three independent
elastic constants of bulk Germanium (diamond cubic structure): **C₁₁**, **C₁₂**, and **C₄₄**.

The method uses **finite differences** (central difference approximation). Each
elastic constant is obtained by applying a small strain ε to the equilibrium cell,
computing the stress tensor σ for both +ε and −ε deformations, and extracting the
elastic constant as:

```
C = [σ(+ε) − σ(−ε)] / (2ε)
```

Both deformations use a strain amplitude of **ε = 0.01** (1%).

---

## Folder Structure

```
week6/
├── Ge.upf               # Germanium pseudopotential (UPF format)
├── C1111/               # Calculation of C₁₁ (uniaxial strain along x)
│   ├── plus_e.in        # Input: ε₁₁ = +0.01
│   ├── minus_e.in       # Input: ε₁₁ = −0.01
│   └── results/
│       ├── plus_e.out   # Output from plus_e.in
│       └── minus_e.out  # Output from minus_e.in
└── C1212/               # Calculation of C₄₄ (shear strain in the 1-2 plane)
    ├── plus_e.in        # Input: ε₁₂ = +0.01
    ├── minus_e.in       # Input: ε₁₂ = −0.01
    └── results/
        ├── plus_e.out   # Output from plus_e.in
        └── minus_e.out  # Output from minus_e.in
```

---

## The Deformations

Both calculations start from the equilibrium diamond-cubic cell of Ge with
lattice parameter **a = 5.76 Å**. The cell is deformed by modifying the
`CELL_PARAMETERS` block.

### C1111 — Uniaxial strain (tetragonal deformation)

The unit cell is uniformly strained along the **x**-direction:

| File        | CELL_PARAMETERS (in units of alat)         |
|-------------|---------------------------------------------|
| `plus_e.in` | `1.01  0.00  0.00 / 0.00  1.00  0.00 / 0.00  0.00  1.00` |
| `minus_e.in`| `0.99  0.00  0.00 / 0.00  1.00  0.00 / 0.00  0.00  1.00` |

This single set of calculations yields **two** elastic constants.

**C₁₁** comes from the longitudinal stress response:

```
C₁₁ = [σ₁₁(+ε) − σ₁₁(−ε)] / (2ε)
```

where σ₁₁ is the (1,1) diagonal component of the stress tensor.

**C₁₂** comes from the transverse stress response to the same deformation:

```
C₁₂ = [σ₂₂(+ε) − σ₂₂(−ε)] / (2ε)
```

where σ₂₂ is the (2,2) diagonal component (equivalently, σ₃₃ by symmetry).

### C1212 — Shear strain (monoclinic deformation)

A shear deformation is applied in the **1-2 plane**:

| File        | CELL_PARAMETERS (in units of alat)          |
|-------------|----------------------------------------------|
| `plus_e.in` | `1.00  0.01  0.00 / 0.00  1.00  0.00 / 0.00  0.00  1.00` |
| `minus_e.in`| `1.00 -0.01  0.00 / 0.00  1.00  0.00 / 0.00  0.00  1.00` |

The elastic constant **C₄₄** is then:

```
C₄₄ = [σ₁₂(+ε) − σ₁₂(−ε)] / (2ε)
```

where σ₁₂ is the (1,2) off-diagonal component of the stress tensor.
For a cubic crystal, C₄₄ = C₅₅ = C₆₆.

---

## Running the Calculations

### Prerequisites

- Quantum ESPRESSO ≥ 6.x installed (the `pw.x` executable must be in your `$PATH`)
- The pseudopotential file `Ge.upf` must be accessible from each input directory

### Step 1 — Copy the pseudopotential

The `pseudo_dir` in each input file is set to `'./'`, so copy `Ge.upf` into the
folder you are running from:

```bash
cp ../../Ge.upf .
```

### Step 2 — Run pw.x

Run both jobs for each elastic constant in sequence (or submit them to the cluster;
see below):

```bash
# Example: run C1111 calculations locally (serial)
cd C1111
cp ../Ge.upf .
pw.x < plus_e.in  > results/plus_e.out
pw.x < minus_e.in > results/minus_e.out
```

Repeat the same steps inside the `C1212/` folder.

### Step 3 — Extract the stress tensor

After the calculation completes, search for the stress block in the output file:

```bash
grep -A 4 "total   stress" results/plus_e.out
```

The output looks like:

```
     total   stress  (Ry/bohr**3)                   (kbar)     P=   -3.32
  -0.00005119  -0.00000000   0.00000000    -7.53  -0.00   0.00
  -0.00000000  -0.00000830   0.00000000    -0.00  -1.22   0.00
   0.00000000   0.00000000  -0.00000830     0.00   0.00  -1.22
```

The right-hand 3×3 block contains the stress tensor components in **kbar**.
- For **C₁₁**: read σ₁₁ (row 1, col 1), e.g. −7.53 kbar
- For **C₁₂**: read σ₂₂ (row 2, col 2), e.g. −1.22 kbar
- For **C₄₄**: read σ₁₂ (row 1, col 2), e.g. −6.79 kbar

### Step 4 — Compute the elastic constants

Apply the central difference formula (ε = 0.01):

```
C = [σ(+ε) − σ(−ε)] / (2 × 0.01)
```

The result will be in **kbar**. Convert to GPa by dividing by 10.

> **Sign convention**: In QE's output, positive diagonal stress components
> correspond to compressive stress, and negative to tensile stress. Keep track
> of signs when computing the difference.

---

## Submitting to the Cluster (SLURM)

Below is a template submission script. Adjust the account, partition, and module
names to match your cluster's configuration.

```bash
#!/bin/bash
#SBATCH --job-name=Ge_elastic
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=<YOUR_PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=01:00:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

module load quantum-espresso   # adjust to your cluster's module name

cp ../../Ge.upf .

mpirun -np $SLURM_NTASKS pw.x < plus_e.in  > results/plus_e.out
mpirun -np $SLURM_NTASKS pw.x < minus_e.in > results/minus_e.out
```

Save this script (e.g., `submit.sh`) inside each deformation folder and submit with:

```bash
sbatch submit.sh
```

---

## Key Parameters

| Parameter   | Value  | Description                              |
|-------------|--------|------------------------------------------|
| `ecutwfc`   | 80 Ry  | Plane-wave kinetic energy cutoff         |
| `degauss`   | 0.01 Ry| Gaussian smearing width                  |
| k-points    | 4×4×4  | Monkhorst–Pack grid (no offset)          |
| `nosym`     | `.true.`| Symmetry disabled (required for finite-difference deformations) |
| ε           | ±0.01  | Strain amplitude (1%)                    |

---

## Expected Results

Pre-computed outputs are available in each `results/` subdirectory and can be used
to verify your analysis workflow before running new jobs.

| Constant | This calculation | Experiment |
|----------|-----------------|------------|
| C₁₁      | ~102 GPa        | 126 GPa    |
| C₁₂      | ~39 GPa         | 44 GPa     |
| C₄₄      | ~68 GPa         | 66 GPa     |

Discrepancies with experiment reflect the choice of pseudopotential, exchange-
correlation functional, and convergence parameters.
