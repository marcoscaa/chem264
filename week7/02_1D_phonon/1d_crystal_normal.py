import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

# --- Physical Parameters ---
a = 1.0           
N = 40            
n = np.arange(N)  
x_eq = n * a      

# Define allowed wavevectors (1st Brillouin Zone)
m_values = np.arange(-N//2 + 1, N//2 + 1)
qs_valid = (2 * np.pi * m_values) / (N * a)

# --- Figure Setup ---
fig, ax = plt.subplots(figsize=(10, 4))
plt.subplots_adjust(bottom=0.25)
fig.canvas.manager.set_window_title('Single Normal Mode')

line, = ax.plot([], [], 'bo', markersize=8, markeredgecolor='black')
ax.set_xlim(-1.5, N * a + 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_yticks([])
ax.set_title("Single Normal Mode $u_n(t) = u \cdot \cos(qna - \omega t)$")

# --- Slider ---
ax_q = plt.axes([0.15, 0.1, 0.7, 0.03])
slider_q = Slider(ax_q, 'Wavevector $q$', qs_valid[0], qs_valid[-1], valinit=qs_valid[N//4], valstep=qs_valid[1]-qs_valid[0])

# --- Animation Logic ---
def animate(frame):
    # SPEEDUP: Time step doubled to 0.30 for 2X animation speed
    t = frame * 0.30      
    q = slider_q.val
    omega = np.abs(np.sin(q * a / 2))
    
    displacements = 0.2 * np.cos(q * x_eq - omega * t)
    line.set_data(x_eq + displacements, np.zeros_like(x_eq))
    return line,

ani = FuncAnimation(fig, animate, frames=200, interval=40, blit=True)

plt.show()