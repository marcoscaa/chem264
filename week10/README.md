# Week 10 — Dielectric function of Ge and C (diamond structure)

This folder contains the Quantum ESPRESSO input files used to compute the
**frequency-dependent dielectric function** $\epsilon(\omega) = \epsilon_1(\omega) + i\epsilon_2(\omega)$
of germanium and carbon in the diamond structure, using the `pw.x` (plane-wave DFT)
and `epsilon.x` (post-processing) codes.

These calculations are the practical companion to the **Week 10 lecture notes**
(section *"Practical calculation of the dielectric response"*). The lecture derives
the dielectric function from time-dependent perturbation theory and Fermi's golden
rule; here we evaluate it numerically from the Kohn–Sham band structure.

---

## 1. Background: what we are computing

In the lecture we showed that, within first-order (linear) response theory, an
external oscillating electric field polarizes the electron density and the response
is governed by the complex dielectric function. Expressed in terms of Kohn–Sham
states it reads (lecture, *"Dielectric susceptibility from KS states"*):

$$
\epsilon_2(\omega) = \frac{4\pi}{V}\int_{\mathrm{BZ}}\frac{d^3k}{(2\pi)^3}
\sum_{m}\sum_{n\neq m} f_m(\mathbf{k})\,|\boldsymbol{\mu}_{nm}(\mathbf{k})|^2
\left(\frac{\Gamma}{(\omega_{nm}(\mathbf{k})-\omega)^2+\Gamma^2} -
\frac{\Gamma}{(\omega_{nm}(\mathbf{k})+\omega)^2+\Gamma^2}\right)
$$

with an analogous expression for $\epsilon_1$. Here the indices $m$ and $n$ run
over **energy states** (Kohn–Sham bands): the occupation number $f_m$ guarantees
that we only count transitions out of *occupied* states $m$ into states $n$
(the occupation depends on $\mathbf{k}$ as well, hence $f_m(\mathbf{k})$). The
$\int_{\mathrm{BZ}} d^3k$ is a genuine integral over reciprocal space, and the
transition energy $\omega_{nm}(\mathbf{k})=E_n(\mathbf{k})-E_m(\mathbf{k})$ depends
on the crystal momentum $\mathbf{k}$. In practice `epsilon.x` evaluates this
Brillouin-zone integral **numerically**, as a weighted sum over the discrete
k-points of the NSCF mesh — the denser the grid, the closer the numerical sum
approaches the true integral.

The key physical ingredients map directly onto the inputs you will run:

| Quantity in the theory | Where it comes from in the calculation |
|---|---|
| KS eigenvalues $E_n(\mathbf{k})$, transition energies $\omega_{nm}(\mathbf{k})=E_n(\mathbf{k})-E_m(\mathbf{k})$ | the SCF + NSCF band structure (`00_scf.in`, `01_nscf.in`) |
| Occupations $f_m(\mathbf{k})$ (which states are valence vs. conduction) | smearing in the SCF/NSCF runs |
| Dipole matrix elements $\|\boldsymbol{\mu}_{nm}(\mathbf{k})\|^2$ | computed internally by `epsilon.x` |
| Lifetime broadening $\Gamma$ | `intersmear` in `02_epsilon.in` |
| Integral over the Brillouin zone $\int_{\mathrm{BZ}} d^3k$ (done numerically as a sum over the k-grid) | the dense **20×20×20** k-mesh in the NSCF run |
| Unit-cell volume $V$ | the lattice from `ibrav`/`a` |

Recall the take-home physics from the lecture: only the **imaginary part**
$\epsilon_2$ corresponds to energy dissipation (light *absorption*), since the
time-averaged power loss is $W_\text{ave}=\tfrac{1}{2}F_0^2\,\omega\,\chi_2(\omega)$.
The **real part** $\epsilon_1$ is the in-phase (non-absorbing) response, and
$\epsilon_1(\omega\to 0)$ gives the static dielectric constant.

---

## 2. The three-step workflow

The dielectric function is obtained with three input files that must be run **in
order**, because each step reuses the output of the previous one (stored under
`outdir`, e.g. `./Ge_out/`).

### Step 1 — `00_scf.in` (self-consistent field)

A standard SCF calculation that converges the Kohn–Sham potential and the ground-state
electron density on a moderate **8×8×8** k-mesh.

- `ibrav = 2`, `a = 5.76` Å (Ge) / `3.56` Å (C): FCC primitive cell of the diamond
  structure, 2 atoms per cell.
- `ecutwfc = 60` Ry: plane-wave cutoff.
- Gaussian `smearing` with a small `degauss = 0.01` Ry: appropriate for
  semiconductors and helps SCF convergence.

### Step 2 — `01_nscf.in` (non-self-consistent field)

Reuses the converged potential from Step 1 and recomputes the bands on a **much
denser 20×20×20 k-mesh**, adding empty (conduction) states. This dense sampling is
needed to accurately perform the Brillouin-zone integral $\int_{\mathrm{BZ}} d^3k$
in the formulas above (which the code evaluates as a numerical sum over this grid). Two changes relative to the SCF run are essential:

- `nbnd = 20`: adds unoccupied bands. Absorption is a transition from occupied to
  empty states, so we need conduction bands as final states $n$.
- `nosym = .true.` and `noinv = .true.`: **disable symmetry**. By default QE uses
  crystal symmetry to evaluate only the inequivalent k-points, but `epsilon.x` does
  **not** support symmetry-reduced meshes — it needs the full, explicit grid.

### Step 3 — `02_epsilon.in` (post-processing with `epsilon.x`)

Reads the KS eigenvalues, occupations and wavefunctions from Step 2 and assembles
$\epsilon_1(\omega)$ and $\epsilon_2(\omega)$. Key parameters in `&energy_grid`:

- `calculation = 'eps'`: compute the dielectric function.
- `wmin`, `wmax`, `nw`: frequency window (0–30 eV) sampled on 600 points.
- `intersmear = 0.15` eV: the broadening $\Gamma$ (interband). Larger values give a
  smoother spectrum but wash out fine structure.
- `intrasmear = 0.0` eV: intraband (Drude) broadening, relevant only for metals.
- `shift`: a rigid shift of the spectrum (e.g. a crude "scissor" correction for the
  DFT band-gap error); left at 0 here.

> **Important:** `prefix` and `outdir` must be **identical** across all three files
> (`Ge_scf` / `./Ge_out/` for germanium, `C_scf` / `./C_out/` for carbon), otherwise
> later steps cannot find the files written by earlier ones. The `outdir` must exist
> before you start (`epsilon.x` will not create it).

---

## 3. How to run

### Option A — interactively (small jobs / local machine)

From inside `Ge_diamond/` (or `C_diamond/`), with Quantum ESPRESSO in your `PATH`:

```bash
mkdir -p Ge_out                          # create outdir if it does not exist
pw.x      < 00_scf.in     > 00_scf.out    # Step 1: SCF
pw.x      < 01_nscf.in    > 01_nscf.out   # Step 2: dense NSCF + empty bands
epsilon.x < 02_epsilon.in > 02_epsilon.out # Step 3: dielectric function
```

(For carbon, replace `Ge_out` with `C_out`.)

### Option B — on the cluster with SLURM

A ready-to-use submission script `sub.cmd` is provided in each folder. Edit the
`--mail-user=<cruzid>@ucsc.edu` line to your CruzID, then submit:

```bash
sbatch sub.cmd
```

The script loads `quantumespresso/7.2`, runs `pw.x` in parallel with k-point
pooling (`-nk 4`), and finishes with the serial `epsilon.x` post-processing — the
same three steps as above.

---

## 4. Output files

`epsilon.x` writes three text files (the `Ge`/`C` in the names follows the `prefix`):

| File | Content |
|---|---|
| `epsi_*.dat` | **imaginary** part $\epsilon_2(\omega)$ — the absorption spectrum |
| `epsr_*.dat` | **real** part $\epsilon_1(\omega)$ |
| `eels_*.dat` | electron energy-loss spectrum (EELS), $-\mathrm{Im}[1/\epsilon]$ |
| `ieps_*.dat` | the inverse dielectric function |

Each file has **4 columns**: column 1 is the energy in **eV**, and columns 2–4 are
the response along the *x*, *y* and *z* directions. For the cubic diamond structure
the three directions are equivalent, so the spectrum is isotropic (the three columns
coincide).

Reference outputs from a completed run are stored in each `results/` subfolder.

### Quick plot

```bash
# Imaginary part (absorption) vs. energy, x-component
gnuplot -p -e "plot 'results/epsi_Ge_scf.dat' u 1:2 w l title 'Im eps (Ge)'"
```

or with Python:

```python
import numpy as np, matplotlib.pyplot as plt
e, x, y, z = np.loadtxt('results/epsi_Ge_scf.dat', unpack=True)
plt.plot(e, x); plt.xlabel('Energy (eV)'); plt.ylabel(r'$\epsilon_2$'); plt.show()
```

---

## 5. What to look for in the results

This is where the calculation connects back to the physics of the lecture:

- **Absorption onset.** The energy where $\epsilon_2$ first rises from zero is set
  by the band gap (the smallest $\omega_{nm}(\mathbf{k})$ with nonzero dipole). As discussed in
  the lecture (Fig. of the dielectric functions), **PBE underestimates band gaps**:
  it predicts **Ge to be nearly metallic** (vanishing gap, $\epsilon_2$ rising almost
  from zero energy), while **diamond C is correctly insulating** with an onset
  around ~4 eV. Comparing the two materials makes the DFT gap problem concrete.
- **Peaks in $\epsilon_2$** correspond to regions of the Brillouin zone where many
  valence→conduction transitions are resonant with $\omega$ (van Hove singularities
  in the joint density of states), weighted by the dipole matrix elements.
- **Static dielectric constant.** Read $\epsilon_1$ at $\omega\to 0$; Ge (more
  polarizable, smaller gap) has a much larger value than diamond C.
- **Effect of broadening.** Re-run Step 3 with a larger `intersmear` and watch the
  peaks smear out — this is the phenomenological lifetime $\Gamma$ from the lecture.

---

## 6. Things to try (exercises)

1. **k-point convergence.** Increase the NSCF mesh (e.g. 24×24×24) and check whether
   the spectrum changes — the BZ integral must be converged.
2. **Number of bands.** Increase `nbnd` and confirm the high-energy part of the
   spectrum (up to `wmax`) is converged with respect to the number of conduction
   states.
3. **Broadening.** Vary `intersmear` (0.05 → 0.3 eV) and relate the result to the
   lifetime broadening $\Gamma$ in the lecture formulas.
4. **Scissor shift.** Use `shift` to rigidly move the spectrum and mimic a corrected
   band gap, then compare with experiment.

---

## 7. Folder contents

```
week10/
├── Ge_diamond/            # Germanium, diamond structure
│   ├── 00_scf.in          # Step 1: SCF
│   ├── 01_nscf.in         # Step 2: dense NSCF + empty bands, symmetry off
│   ├── 02_epsilon.in      # Step 3: epsilon.x post-processing
│   ├── Ge.upf             # pseudopotential
│   ├── sub.cmd            # SLURM submission script
│   └── results/           # reference outputs (.out and .dat files)
└── C_diamond/             # Carbon (diamond), same workflow
    ├── 00_scf.in
    ├── 01_nscf.in
    ├── 02_epsilon.in
    ├── C.upf
    ├── sub.cmd
    └── results/
```
