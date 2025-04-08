import numpy as np
import scipy.linalg as lina
import plot
from potentials import softened_coulomb, softened_mean_field_repulsion
from utils import make_pos, discretize_space

def H2_mean_field(pos, x, N, electron_density, independent=False):

    dx = x[1] - x[0]
    num_points = x.size

    # Potential energy (ion-electron attraction)
    V = 0
    for pos_ions in pos:
        V += softened_coulomb(x, pos_ions)

    if (not independent):
        #Potential energy (electron-electron mean field repulsion)
        V += softened_mean_field_repulsion(electron_density, x, dx) 

    # Kinetic energy operator (finite difference)
    T = -0.5 * (-2 * np.eye(num_points) + np.eye(num_points, k=1) + np.eye(num_points, k=-1)) / dx**2

    # Hamiltonian
    H = T + np.diag(V)

    # Solve for eigenvalues and eigenvectors
    eigenvalues, eigenvectors = lina.eigh(H)

    return eigenvalues, eigenvectors, x

def compute_electron_density(eigenvectors,N):
    density = 2*np.sum(np.abs(eigenvectors[:,:int(N/2)])**2,axis=1)
    return density

def scf(pos, x, N, maxiter=10):

    #Initial guess using the independent electron approximation
    eigenvalues, eigenvectors, x = H2_mean_field(pos, x, N, None, independent=True) 

    for it in range(maxiter):
        electron_density = compute_electron_density(eigenvectors,N) 
        eigenvalues, eigenvectors, x = H2_mean_field(pos, x, N, electron_density) 
        energy = compute_energy(eigenvalues,electron_density,N,pos,x)
        print(it, energy)

    return eigenvalues, eigenvectors, x

def compute_energy(eigenvalues,electron_density,N,pos,x):
    V_rep = softened_mean_field_repulsion(electron_density, x, x[1]-x[0])
    E_rep = np.trapz(electron_density*V_rep,x)
    E_nn = 0.0
    for I in range(N):
        for J in range(I,N):
            E_nn -= softened_coulomb(pos[I],pos[J])
    E = 2*np.sum(eigenvalues[:N//2]) - E_rep + E_nn 
    return E

def main():

    N=2
    R = 1.4  # Internuclear distance (Bohr radii)
    pos = make_pos(N,R)
    x = discretize_space(pos)

    #Solve equations self-consistently
    eigenvalues, eigenvectors, x = scf(pos, x, N, 4)
    energy = compute_energy(eigenvalues,eigenvectors,N,pos,x)

    plot.plot_potential(pos,x)
    plot.plot_eigenvectors(x,eigenvectors,eigenvalues,N)
    print(eigenvalues[:N])

main()
