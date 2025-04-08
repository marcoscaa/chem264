import numpy as np
import scipy.linalg as lina
import plot
from potentials import softened_coulomb
from utils import make_pos, discretize_space

def H2_independent_electron(pos, x, N):

    dx = x[1] - x[0]
    num_points=x.size

    # Potential energy (independent electron approximation)
    V = 0
    for pos_n in pos:
        V += softened_coulomb(x, pos_n)

    # Kinetic energy operator (finite difference)
    T = -0.5 * (-2 * np.eye(num_points) + np.eye(num_points, k=1) + np.eye(num_points, k=-1)) / dx**2

    # Hamiltonian
    H = T + np.diag(V)

    # Solve for eigenvalues and eigenvectors
    eigenvalues, eigenvectors = lina.eigh(H)

    return eigenvalues, eigenvectors, x

def main():

    N=4
    R = 1.4  # Internuclear distance (Bohr radii)
    pos = make_pos(N,R) 
    x = discretize_space(pos,0.01) 

    eigenvalues, eigenvectors, x = H2_independent_electron(pos, x ,N)

    plot.plot_potential(pos,x)
    plot.plot_eigenvectors(x,eigenvectors,eigenvalues,N)
    print(eigenvalues[:N])

#main()
