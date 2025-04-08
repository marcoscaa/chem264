import numpy as np
import matplotlib.pyplot as plt
from potentials import softened_coulomb

def plot_potential(pos,x):
    V=0
    for pos_ions in pos:
        Vp = softened_coulomb(x, pos_ions)
        plt.plot(x, Vp, label='Nucleus {n+1} Potential')
        V+=Vp
    plt.plot(x, V, label='Combined Potential')
    plt.xlabel('x (Bohr)')
    plt.ylabel('Potential Energy (Hartree)')
    plt.xlim([np.min(x)+2,np.max(x)-2])
    plt.legend()
    plt.show()

def plot_eigenvectors(x,eigenvectors,eigenvalues,N):
    for n in range(N):
        plt.plot(x, eigenvectors[:, n], label=f'MO {n+1} (E = {eigenvalues[n]:.3f})')
    plt.xlabel('x (Bohr radii)')
    plt.ylabel('Wave Function')
    plt.xlim([np.min(x)+2,np.max(x)-2])
    plt.legend()
    plt.show()

def plot_eigenvalues(eigenvalues):
    plt.plot(eigenvalues)
    plt.ylabel("Eigenvalue (Hartree)")
    plt.xlabel("State index")
    plt.show()

