# Week 3 — First-Principles Molecular Dynamics of a Water Molecule

This activity introduces Born-Oppenheimer Molecular Dynamics (BOMD) using [Quantum ESPRESSO](https://www.quantum-espresso.org/). You will simulate a single H₂O molecule in a cubic box at constant energy (NVE ensemble) and observe how the atomic positions and velocities evolve under the forces computed from DFT at each time step.

## System

**Molecule:** Single H₂O (3 atoms) in a 10 Å cubic supercell (explicit `CELL_PARAMETERS`)  
**Pseudopotentials:** Norm-conserving, scalar-relativistic, PBE (`H.upf`, `O.upf` — PseudoDojo ONCVPSP)  
**Functional:** PBE (Perdew-Burke-Ernzerhof, GGA)  
**Ensemble:** NVE (constant number of particles, volume, and energy — no thermostat)  
**Time step:** 10 Rydberg a.u. ≈ 0.48 fs  
**Total simulation time:** 2000 steps ≈ ~1 ps

---

## Getting Started on the Cluster

### 1. Log in

Connect to the UCSC HummingBird cluster with X11 forwarding enabled (needed for graphical tools):

```bash
ssh -XY <cruzid>@hb.ucsc.edu
```

Replace `<cruzid>` with your UCSC CruzID.

### 2. Get the course repository

If you have already cloned the repository in a previous week, update it:

```bash
cd ~/chem264
module load git
git pull
```

If this is your first time, clone it first:

```bash
module load git
git clone https://github.com/marcoscaa/chem264.git ~/chem264
cd ~/chem264
```

### 3. Navigate to this week's folder

```bash
cd ~/chem264/codes/week3/FPMD
```

---

## Files in This Directory

```
week3/FPMD/
├── h2o_nvt.in    # pw.x input file for BOMD (NVE)
├── sub.cmd       # SLURM job submission script
└── pseudo/       # Directory where pseudopotential files must be placed
    ├── O.upf
    └── H.upf
```

The pseudopotential files (`O.upf` and `H.upf`) should already be present in the `pseudo/` directory. They were downloaded from [pseudo-dojo.org](https://www.pseudo-dojo.org) (NC SR PBE standard table, UPF2 format).

---

## Visualizing the Input File

Inspect the pw.x input file before submitting:

```bash
cat h2o_nvt.in
```

Or open it in a text editor:

```bash
vi h2o_nvt.in
```

We can also look at the atomic structure of the system using ase:

```bash
module load minoconda3
conda activate chem264
ase gui h2o_nvt.in 
```

### Input File Walkthrough

```fortran
&CONTROL
  calculation   = 'md'          ! Run BOMD (molecular dynamics)
  restart_mode  = 'from_scratch'
  pseudo_dir    = './pseudo'    ! Directory containing pseudopotential files
  outdir        = './tmp'       ! Temporary files go here (wavefunctions, charge density)
  prefix        = 'h2o'         ! Prefix for all output files
  dt            = 10.0          ! Time step in Rydberg a.u. (~0.48 fs)
  nstep         = 2000          ! Number of MD steps (~1 ps total)
  iprint        = 1             ! Write energies/forces at every MD step
/
```

The main difference from the SCF runs you did in week 2 is `calculation = 'md'`. This tells pw.x to perform molecular dynamics instead of a static self-consistent field calculation. At each MD step, pw.x solves the Kohn-Sham equations to get the ground-state electron density, computes the forces on the nuclei, and then moves the atoms according to Newton's equations of motion. `dt` sets the nuclear time step and `nstep` the total number of steps.

```fortran
&SYSTEM
  ibrav         = 0             ! Cell vectors given explicitly in CELL_PARAMETERS
  nat           = 3             ! 3 atoms: 1 O + 2 H
  ntyp          = 2             ! 2 atomic species: O and H
  ecutwfc       = 60.0          ! Plane-wave kinetic energy cutoff (Ry)
  nosym         = .true.        ! Disable symmetry — required for MD
/
```

This section is largely the same as in the SCF case. Setting `ibrav = 0` tells pw.x that the cell vectors are provided explicitly in the `CELL_PARAMETERS` card below, rather than being generated from a Bravais-lattice index and lattice parameters. `nosym = .true.` is important for MD: as atoms move, the symmetry of the system changes at every step, so symmetry analysis must be turned off to avoid errors.

Note that `ecutrho` is not set here — for norm-conserving pseudopotentials, pw.x defaults to `ecutrho = 4 * ecutwfc`, which is the correct value.

```fortran
&ELECTRONS
  electron_maxstep = 100
  conv_thr         = 1.0d-8    ! Tight convergence to ensure accurate forces
/
```

The electrons section is identical to what you used in the SCF runs. A tight convergence threshold is important in MD because the nuclear forces are the gradient of the total energy — if the electron density is not well converged, the forces will be noisy and the trajectory unreliable.

```fortran
&IONS
  ion_dynamics      = 'verlet'          ! Velocity-Verlet integrator (NVE)
  ion_velocities    = 'from_input'      ! Velocities read from VELOCITY section
  ion_temperature   = 'not_controlled'  ! No thermostat — NVE ensemble
  pot_extrapolation = 'second_order'    ! Extrapolate V_eff from previous steps
  wfc_extrapolation = 'second_order'    ! Extrapolate wavefunctions from previous steps
/
```

This section is new and controls the nuclear dynamics:

- `ion_dynamics = 'verlet'` selects the velocity-Verlet algorithm, which propagates atomic positions and velocities and exactly conserves the total energy (NVE).
- `ion_velocities = 'from_input'` tells pw.x to read the initial velocities from the `ATOMIC_VELOCITIES` card in the input file, rather than setting them to zero or randomizing them internally.
- `ion_temperature = 'not_controlled'` means no thermostat is applied — the system evolves at constant total energy.
- `pot_extrapolation` and `wfc_extrapolation` are efficiency settings: instead of starting the SCF from scratch at every MD step, pw.x extrapolates the potential and wavefunctions from the previous steps as an initial guess, reducing the number of SCF iterations needed.

Note: in QE pw.x 7.2, the `tempw` keyword only initializes velocities when an active thermostat is chosen. With `ion_temperature = 'not_controlled'` it is ignored and atoms start at rest. Initial velocities must therefore be supplied explicitly via the `ATOMIC_VELOCITIES` card, and `ion_velocities = 'from_input'` must be set to instruct pw.x to read them.

```fortran
ATOMIC_SPECIES
  O   15.9994  O.upf
  H    1.0079  H.upf

CELL_PARAMETERS {angstrom}
  10.0   0.0   0.0
   0.0  10.0   0.0
   0.0   0.0  10.0

ATOMIC_POSITIONS {angstrom}
  O   5.000000   5.000000   5.000000
  H   5.757800   5.585100   5.000000
  H   4.242200   5.585100   5.000000

ATOMIC_VELOCITIES
  O  -0.0002572821  -0.0000543519   0.0000927743
  H   0.0019967664  -0.0003758507  -0.0005466765
  H   0.0020873291   0.0012386328  -0.0009260225

K_POINTS {gamma}
```

The `CELL_PARAMETERS` card explicitly defines the three lattice vectors of the simulation box in Å, forming a 10 Å cubic supercell. This replaces the `ibrav`/`celldm` approach and makes the cell geometry immediately visible in the input file. The molecule is placed at the center of the box using the experimental geometry (O-H bond = 0.9572 Å, H-O-H angle = 104.52°). Since this is a gas-phase molecule (non-periodic), only the Γ-point is needed for k-point sampling.

The `ATOMIC_VELOCITIES` card provides the initial nuclear velocities in Rydberg atomic units (Bohr/(ħ/Ry)). This is the time unit used internally by pw.x: `dt = 10` corresponds to 10 ħ/Ry ≈ 0.48 fs, confirming the unit. Each velocity component is drawn from a Gaussian with zero mean and standard deviation σ = √(2 k_B T / m), which follows from the equipartition theorem (½ k_B T per degree of freedom) expressed in Rydberg a.u., where the kinetic energy is ¼mv² rather than the familiar ½mv². After sampling, the center-of-mass velocity is subtracted to remove net translation, and the velocities are uniformly rescaled so that the instantaneous temperature — defined as KE / (½ N_dof k_B), with N_dof = 3N − 3 = 6 — is exactly 300 K.

---

## Visualizing the Submission Script

```bash
cat sub.cmd
```

```bash
#!/bin/bash

#SBATCH -p instruction          # Submit to the instructional partition
#SBATCH -J h2o_md               # Job name shown in the queue
#SBATCH --mail-user=<cruzid>@ucsc.edu
#SBATCH --mail-type=ALL         # Email on job start, end, and failure
#SBATCH -o job%.j.out           # Redirect stdout to this file
#SBATCH -n 4                    # Request 4 MPI tasks
#SBATCH -t 01:00:00             # Wall-clock time limit: 1 hour
#SBATCH --mem=1G                # Memory per node

export OMPI_MCA_btl=tcp,sm,self
module load quantumespresso/7.2

mpirun -np $SLURM_NTASKS pw.x < h2o_nvt.in > h2o_nvt.out
```

The script uses SLURM directives (`#SBATCH`) to request resources on the cluster. The key lines are:

- `-p instruction` — submits to the instructional queue reserved for this course.
- `-n 4` — requests 4 MPI processes; pw.x will parallelize the calculation across them.
- `-t 01:00:00` — requests a 1-hour time limit. The ~1 ps BOMD run should complete well within this.
- `module load quantumespresso/7.2` — loads the Quantum ESPRESSO module on the cluster.
- The last line runs pw.x with MPI, reading from `h2o_nvt.in` and writing all output to `h2o_nvt.out`.

---

## Submitting the Job

Before submitting, edit `sub.cmd` and replace `<cruzid>` with your actual CruzID:

```bash
nano sub.cmd
```

Then create the temporary output directory and submit:

```bash
mkdir -p tmp
sbatch sub.cmd
```

Monitor the job status with:

```bash
squeue -u <cruzid>
```

Once the job starts, you can follow the output in real time:

```bash
tail -f h2o_nvt.out
```

Each MD step will print the step number, total energy, kinetic energy of the nuclei, and temperature. A well-behaved NVE run should show a conserved total energy (kinetic + potential) throughout the simulation.

---

## Analyzing the Results

Once the job has finished (or while it is still running), you can plot the energies and temperature as a function of time using the provided Python script.

First, activate the course Python environment:

```bash
module load miniconda3
conda activate chem264
```

Then run the analysis script, passing the output file as an argument:

```bash
python plot_md.py h2o_nvt.out
```

The script will open an interactive plot window with two panels:

1. **Energies vs time** — potential energy $E_\text{pot}$, nuclear kinetic energy $E_\text{kin}$, and the conserved total energy $E_\text{tot} = E_\text{pot} + E_\text{kin}$, all in eV. The energies are shifted so the mean total energy is zero, making fluctuations easy to read.
2. **Temperature vs time** — instantaneous nuclear temperature in K, with the time-averaged mean shown as a dashed line.

The figure is also saved to `md_energies.png` in the current directory. In a well-behaved NVE simulation, $E_\text{tot}$ should remain nearly flat throughout, while $E_\text{pot}$ and $E_\text{kin}$ oscillate out of phase as kinetic and potential energy are exchanged during the molecular vibrations.
