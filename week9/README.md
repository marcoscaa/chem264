# Band structure and density of states of Ge (diamond structure)

This folder contains the Quantum ESPRESSO (PW) input files needed to compute the
**electronic band structure** and the **density of states (DOS)** of germanium in
the diamond (FCC, two-atom basis) structure. The calculations follow the workflow
described in the *Week 9 — Band structure* lecture notes (`../../Week 9.tex`); the
sections below point back to the relevant theory as you go.

## Background: what we are computing and why

In the lecture we used **Bloch's theorem**,

$$\psi_{n\mathbf{k}}(\mathbf{r}) = u_{n\mathbf{k}}(\mathbf{r})e^{i\mathbf{k}\cdot\mathbf{r}}$$

to rewrite the Kohn-Sham (KS) equations for a periodic crystal as one eigenvalue
problem *per* reciprocal-space point $\mathbf{k}$:

$$\left[-\tfrac{1}{2}(\nabla + i\mathbf{k})^2 + V_{tot}(\mathbf{r})\right] u_{n\mathbf{k}}(\mathbf{r}) = \epsilon_{n\mathbf{k}} u_{n\mathbf{k}}(\mathbf{r}).$$

The **band structure** is the plot of the eigenvalues $\epsilon_{n\mathbf{k}}$ along a
path through high-symmetry points of the Brillouin zone, and the **DOS** $g(E)$
counts how many of those eigenvalues fall in each energy window. Both are obtained
in two stages, exactly as in the notes (Sec. *Practical calculations of the band
structure*):

1. **One** self-consistent (SCF) calculation to converge the electron density
   $n(\mathbf{r})$ and the total potential $V_{tot}(\mathbf{r})$ on a uniform
   $\mathbf{k}$ mesh.
2. **Non-self-consistent (nscf / bands)** calculations that reuse the converged
   $V_{tot}$ and simply diagonalize the KS Hamiltonian at the $\mathbf{k}$ points
   we care about — a *path* for the bands, and a *dense uniform grid* for the DOS.

## Files in this folder

| File | Code | Purpose |
|------|------|---------|
| `00_scf.in`      | `pw.x`    | SCF calculation → converged $n(\mathbf{r})$ and $V_{tot}$ |
| `01_bands.in`    | `pw.x`    | `bands` (non-SCF) run: eigenvalues along the high-symmetry path |
| `02_pp_bands.in` | `bands.x` | Post-processing: extract eigenvalues into a plottable file |
| `03_dos.in`      | `pw.x`    | `nscf` run on a dense $20\times20\times20$ grid for the DOS |
| `04_pp_dos.in`   | `dos.x`   | Post-processing: build $g(E)$ from the dense-grid eigenvalues |
| `Ge.upf`         | —         | PBE pseudopotential for germanium (UPF format) |
| `plot_bands.py`  | Python    | Plot the band structure from `ge_bands.dat.gnu` |
| `plot_dos.py`    | Python    | Plot the DOS from `ge_dos.dat` |
| `sub.cmd`        | SLURM     | Batch script that runs all five steps in order |
| `results/`       | —         | Reference output of a completed run (for comparison) |

The **`prefix` (`Ge_scf`) and `outdir` (`./Ge_out/`) must be identical in every
input file** — this is how the non-SCF steps find the density and potential saved
to disk by the SCF run.

## Step-by-step instructions

The five calculations must be run **in order**, because each non-SCF step depends
on the converged charge density written by step 1 into `./Ge_out/`. The
pseudopotential `Ge.upf` is already in this folder (`pseudo_dir = './'`), and
`pw.x` creates the output directory `Ge_out/` automatically on the first run, so
there is nothing to set up beforehand.

### 1. SCF calculation — `00_scf.in`

This converges $n(\mathbf{r})$ and $V_{tot}(\mathbf{r})$ on a uniform
$8\times8\times8$ Monkhorst-Pack mesh (the `K_POINTS {automatic}` block). A few
points worth noting, mirroring the lecture:

- `ibrav = 0` with an explicit `CELL_PARAMETERS` block fixes the FCC primitive
  vectors and the two-atom basis (Ge at `0,0,0` and `1/4,1/4,1/4`) using the
  convention of **Setyawan & Curtarolo** — this is what makes the high-symmetry
  $\mathbf{k}$ coordinates in step 2 correct.
- `ecutwfc = 60 Ry` is the plane-wave kinetic-energy cutoff; `a = 5.76 Å` is the
  relaxed lattice parameter.
- A small Gaussian smearing (`degauss = 0.01 Ry`) is used purely to stabilize
  convergence; Ge is a semiconductor, not a metal.

```bash
mpirun -np 4 pw.x -nk 4 < 00_scf.in > 00_scf.out
```

After this finishes, look in `00_scf.out` for the line

```
the Fermi energy is     8.7649 ev
```

This Fermi energy ($\epsilon_F$, the energy of the highest occupied state at
$T=0$; see the notes) is hard-coded into both plotting scripts — **update it there
if you change the cell, cutoff, or pseudopotential.**

### 2. Bands (non-SCF) calculation — `01_bands.in`

This is the `calculation = 'bands'` run. It does **not** redo the SCF cycle; it
reuses the converged potential and only diagonalizes the KS Hamiltonian at the
$\mathbf{k}$ points along the path

$$\Gamma \to X \to U \to L \to \Gamma \to K \to W \to X,$$

requesting 20 interpolated points between each pair (`K_POINTS {crystal_b}`, in
units of the reciprocal primitive lattice). These are the standard FCC
high-symmetry points from Setyawan & Curtarolo.

```bash
mpirun -np 4 pw.x -nk 4 < 01_bands.in > 01_bands.out
```

### 3. Extract the bands — `02_pp_bands.in`

`bands.x` reads the eigenvalues from `./Ge_out/` and writes them in a tidy,
plottable form. It produces `ge_bands.dat`, `ge_bands.dat.gnu` (the file the
Python script reads), and `ge_bands.dat.rap` (symmetry labels).

```bash
bands.x < 02_pp_bands.in > 02_pp_bands.out
```

The cumulative $x$-axis coordinate of each high-symmetry point is printed in
`02_pp_bands.out`; those numbers are the keys in the `hs_pts` dictionary inside
`plot_bands.py`. If you change the path, update that dictionary.

### 4. Dense non-SCF run for the DOS — `03_dos.in`

The DOS needs a *much* finer, uniform sampling of the Brillouin zone than the SCF
mesh — here a $20\times20\times20$ grid (`calculation = 'nscf'`). More $\mathbf{k}$
points means more eigenvalues, and since $g(E)$ is essentially a histogram of
eigenvalues, a denser grid gives a smoother curve.

```bash
mpirun -np 4 pw.x -nk 4 < 03_dos.in > 03_dos.out
```

### 5. Build the DOS — `04_pp_dos.in`

`dos.x` turns the dense-grid eigenvalues into $g(E)$ over the energy window
`emin = -5`, `emax = 20` eV, writing `ge_dos.dat` (columns: energy, DOS,
integrated DOS).

```bash
dos.x < 04_pp_dos.in > 04_pp_dos.out
```

### 6. Plot the results

```bash
python plot_bands.py    # reads ge_bands.dat.gnu
python plot_dos.py      # reads ge_dos.dat
```

Both scripts draw a line at the Fermi energy (8.7649 eV). In the band plot you
should see the valence bands below $\epsilon_F$ and the conduction bands above it;
in the DOS plot the shaded region (states below $\epsilon_F$) integrates to the
number of valence electrons (28 per cell). The small gap between the filled and
empty states is what makes Ge a **semiconductor** (see the notes, Sec. *Metals,
insulators and semiconductors*).

> **A note on the band gap.** The gap you read off these plots is the *Kohn-Sham*
> gap $E_g^{DFT} = \epsilon_{CBM} - \epsilon_{VBM}$, which with PBE
> **underestimates** the experimental (quasiparticle) gap of Ge — PBE even tends
> to close it almost entirely. This is the *band gap problem* discussed at the end
> of the lecture: the exact $E_g^{qp} = E_g^{DFT} + \Delta V_{xc}$ differs by the
> derivative discontinuity $\Delta V_{xc}$, which LDA/GGA functionals miss.
> Reproducing the experimental gap requires hybrid functionals or the GW method.

## Running everything at once on a cluster

The SLURM script `sub.cmd` runs all five steps in sequence. Before submitting:

1. Replace `<cruzid>` in the `--mail-user` line with your CruzID.
2. Adjust `-N`/`-n` (nodes / MPI tasks) and the `-nk` pool size if needed; the
   `pw.x` calls use `-nk 4` to split work over 4 k-point pools, so keep the number
   of MPI tasks a multiple of that.

```bash
sbatch sub.cmd
```

Submit it on the instructional partition (`-p instruction`); the script loads
`quantumespresso/7.2` automatically.

> **Heads-up:** the DOS line in `sub.cmd` invokes `$PW`, which is not defined in
> the script. If that step fails, replace `$PW` with `pw.x`:
> `mpirun -np $SLURM_NTASKS pw.x -nk 4 -in 03_dos.in > 03_dos.out`.

## Expected results

A completed reference run is provided in `results/` (output files plus
`ge_bands.dat*` and `ge_dos.dat`) so you can check your numbers — in particular,
the Fermi energy (8.7649 eV) and the overall shape of the bands and DOS.

## References

- W. Setyawan and S. Curtarolo, *High-throughput electronic band structure
  calculations: Challenges and tools*, Comput. Mater. Sci. **49**, 299 (2010) —
  source of the FCC high-symmetry $\mathbf{k}$ points used in `01_bands.in`.
  <https://www.sciencedirect.com/science/article/pii/S0927025610002697>
- Course lecture notes (Week 9) (*Band structure*).
