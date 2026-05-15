import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

# --- Physical Parameters ---
a = 1.0           
N = 40            
n = np.arange(N)  
x_eq = n * a      

# Define allowed wavevectors and frequencies
m_values = np.arange(-N//2 + 1, N//2 + 1)
qs_valid = (2 * np.pi * m_values) / (N * a)
omegas_all = np.abs(np.sin(qs_valid * a / 2))

# Random phases for thermal equilibrium
np.random.seed(42) 
random_phases = np.random.uniform(0, 2 * np.pi, size=len(qs_valid))

# Amplitude scaling to prevent atoms from crossing (Lindemann criterion)
U_SUPER_SCALE = 0.02 

# --- Figure Setup ---
fig, ax = plt.subplots(figsize=(10, 4))
plt.subplots_adjust(bottom=0.25)
fig.canvas.manager.set_window_title('Thermal Superposition')

line, = ax.plot([], [], 'ro', markersize=8, markeredgecolor='black', alpha=0.8)
ax.set_xlim(-1.5, N * a + 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_yticks([])
ax.set_title("Thermal Equilibrium Superposition (Lindemann Bounded)")

# --- Slider ---
ax_T = plt.axes([0.15, 0.1, 0.7, 0.03])
slider_T = Slider(ax_T, 'Temperature $T$', 0.0, 5.0, valinit=1.0)

# --- Animation Logic ---
def animate(frame):
    # SPEEDUP: Time step doubled to 0.30 for 2X animation speed
    t = frame * 0.30
    T = slider_T.val
    
    if T < 1e-3:
        weights = np.zeros_like(qs_valid)
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            weights = np.sqrt(T / np.maximum(omegas_all, 1e-6))
        weights[omegas_all < 1e-6] = 0 
    
    matrix = np.cos(qs_valid[:, np.newaxis] * x_eq - omegas_all[:, np.newaxis] * t + random_phases[:, np.newaxis])
    displacements = U_SUPER_SCALE * np.dot(weights, matrix) 
    
    line.set_data(x_eq + displacements, np.zeros_like(x_eq))
    return line,

ani = FuncAnimation(fig, animate, frames=200, interval=40, blit=True)

plt.show()