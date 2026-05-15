# 1D phonon visualizations

Two interactive Python demos that illustrate the concepts derived in the
"Lattice vibrations of a 1D crystal" and "Phonons in a 3D crystal lattice"
sections of the Week 7 lecture notes.

Both scripts use a **monoatomic 1D chain** of `N = 40` atoms with lattice
constant `a = 1`, nearest-neighbor harmonic coupling, and periodic boundary
conditions. With $C/M$ rescaled to 1, the textbook dispersion relation

$$
\omega(q) = 2\sqrt{C/M}\,\left|\sin(qa/2)\right|
$$

reduces to $\omega(q) = |\sin(qa/2)|$, which is what both scripts use. The
allowed wavevectors in the first Brillouin zone follow from Born–von Kármán
boundary conditions,

$$
q_m = \frac{2\pi m}{N a}, \qquad m = -\tfrac{N}{2}+1, \ldots, \tfrac{N}{2}
$$

giving exactly `N = 40` independent modes.

---

## `1d_crystal_normal.py` — a single normal mode

This script visualizes **one phonon at a time**: a traveling wave

$$
u_n(t) = u_0 \cos(q\,n a - \omega(q)\, t)
$$

displacing the atoms of the chain. A slider at the bottom lets you pick the
wavevector $q$ from the discrete set of allowed values in the first Brillouin
zone, and the animation updates in real time.

What to look for:

- **Small $|q|$** (long wavelength): atoms move almost coherently with their
  neighbors. The wave clearly travels along the chain — this is an **acoustic
  mode**. The frequency $\omega(q) \approx |q|a/2$ is small (linear sound-like
  branch).
- **$|q| \to \pi/a$** (Brillouin-zone edge): neighboring atoms move in
  opposite directions. $\omega(q)$ reaches its maximum and the group velocity
  $d\omega/dq \to 0$ — the standing-wave limit.
- **$q = 0$**: all atoms move together; this is the rigid-translation mode
  with $\omega = 0$.

---

## `1d_crystal_thermal.py` — thermal superposition of all modes

The same chain, but now we **superpose all $N$ normal modes simultaneously**
with random phases:

$$
u_n(t) = \sum_{q} A_q\,\cos\!\big(q\,n a - \omega(q)\,t + \varphi_q\big)
$$

The random phases $\varphi_q$ mimic the incoherent thermal excitation of the
lattice, and each mode is weighted by an amplitude that grows with temperature
and (more strongly) with the inverse of its frequency:

$$
A_q \;\propto\; \sqrt{T/\omega(q)}
$$

The $q=0$ rigid-translation mode is excluded (it would cause the whole chain
to drift). A small prefactor (`U_SUPER_SCALE = 0.02`) keeps the per-atom RMS
displacement below ~10% of $a$ (the Lindemann melting criterion), so that the
harmonic picture remains visually self-consistent and atoms do not cross.

A slider lets you change the temperature $T$ from 0 to 5 (in arbitrary units).

What to look for:

- **$T = 0$**: all weights vanish, the chain sits at rest.
- **Small $T$**: low-frequency (small-$q$) modes dominate the motion because
  of the $1/\sqrt{\omega(q)}$ weighting. The chain looks like it is undulating
  with long-wavelength waves.
- **Larger $T$**: short-wavelength modes contribute too. The overall RMS
  displacement grows like $\sqrt{T}$, reproducing the classical equipartition
  scaling.

This is the qualitative picture of phonon-mediated thermal motion in a real
crystal, condensed to one spatial dimension for clarity.

---

## How to run

Both scripts are self-contained — no DFT inputs, no output files, only NumPy
and Matplotlib are required.

```bash
pip install numpy matplotlib       # if not already installed
python 1d_crystal_normal.py        # single-mode demo
python 1d_crystal_thermal.py       # thermal-superposition demo
```

Drag the slider in either window to change $q$ (for the normal-mode demo) or
$T$ (for the thermal demo). Close the figure window to exit.

> The animation uses a Matplotlib backend with GUI support (TkAgg, Qt5Agg,
> MacOSX, …). Running over an ssh session into the cluster will not display
> anything unless X-forwarding is enabled (`ssh -Y`). For best results, run
> these scripts on your local machine.
