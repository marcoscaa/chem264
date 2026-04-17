import sys
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
outfile = sys.argv[1] if len(sys.argv) > 1 else "h2o_nvt.out"
dt_au   = 10.0          # time step in Rydberg a.u. (must match h2o_nvt.in)
RY2EV   = 13.6057       # 1 Ry = 13.6057 eV
AU2FS   = 0.04837772    # 1 Rydberg a.u. of time = 0.04837772 fs

# ------------------------------------------------------------
# Parse output file
# ------------------------------------------------------------
epot, ekin, etot, temp = [], [], [], []

with open(outfile) as f:
    for line in f:
        if "!    total energy" in line:
            epot.append(float(line.split()[-2]))          # Ry
        elif "kinetic energy (Ekin)" in line:
            ekin.append(float(line.split()[-2]))          # Ry
        elif "temperature" in line and "K" in line:
            temp.append(float(line.split()[-2]))          # K
        elif "Ekin + Etot (const)" in line:
            etot.append(float(line.split()[-2]))          # Ry

epot = np.array(epot) * RY2EV
ekin = np.array(ekin) * RY2EV
etot = np.array(etot) * RY2EV
temp = np.array(temp)

nsteps = len(etot)
time   = np.arange(1, nsteps + 1) * dt_au * AU2FS / 1000  # ps

# Shift energies so that the mean total energy is zero (easier to read)
epot -= np.mean(etot)
etot -= np.mean(etot)

# ------------------------------------------------------------
# Plot 1: energies vs time
# ------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(5, 5), sharex=True)

ax = axes[0]
ax.plot(time, epot, label=r'$E_\mathrm{pot}$', linewidth=0.8, color='steelblue')
ax.plot(time, ekin, label=r'$E_\mathrm{kin}$', linewidth=0.8, color='tomato')
ax.plot(time, etot, label=r'$E_\mathrm{tot}$ (conserved)', linewidth=1.0,
        color='black', linestyle='--')
ax.set_ylabel('Energy (eV)')
ax.legend(fontsize=8, frameon=False)
ax.axhline(0, linewidth=0.5, color='gray', linestyle=':')

# ------------------------------------------------------------
# Plot 2: temperature vs time
# ------------------------------------------------------------
ax = axes[1]
ax.plot(time, temp, linewidth=0.8, color='darkorange')
ax.axhline(np.mean(temp), linewidth=1.0, color='black', linestyle='--',
           label=f'mean = {np.mean(temp):.0f} K')
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Temperature (K)')
ax.legend(fontsize=8, frameon=False)

plt.tight_layout()
plt.savefig('md_energies.png', dpi=150)
plt.show()

print(f"Steps parsed       : {nsteps}")
print(f"Total time         : {time[-1]:.3f} ps")
print(f"Mean temperature   : {np.mean(temp):.1f} K")
print(f"E_tot std dev      : {np.std(etot)*1000:.2f} meV  (energy conservation)")
