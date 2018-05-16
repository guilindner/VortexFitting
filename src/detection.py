#!/usr/bin/env/ python
"""Different methods for vortex detection
"""
import numpy as np

import tools

def calc_swirling(a):
    """
    2D Swirling strength
    """
    print("Detection method: 2D swirling strength")
    A = np.zeros((a.u.size, 3, 3))
    A = np.array(
        [[a.derivative['dudx'].ravel(), a.derivative['dudy'].ravel(),
          a.derivative['dudz'].ravel()],
         [a.derivative['dvdx'].ravel(), a.derivative['dvdy'].ravel(),
          a.derivative['dvdz'].ravel()],
         [a.derivative['dwdx'].ravel(), a.derivative['dwdy'].ravel(),
          -a.derivative['dudx'].ravel()-a.derivative['dvdy'].ravel()]])

    A = A.transpose(2, 1, 0)
    eigenvalues = np.linalg.eigvals(A)
    swirling = np.max(eigenvalues.imag, axis=1).reshape(a.u[:, 0].size, a.u[0, :].size)
    print('Max value of swirling:', np.max(swirling))
    return swirling

def q_criterion(a):
    """
    Q Criterion
    vorticity magnitude and mean strain rate
    """
    print("Detection method: Q criterion")
    Q = np.zeros((a.u.shape[0], a.u.shape[1]))
    print(a.u.shape[0], a.u.shape[1])
    print(Q.shape)
    for i in range(a.u.shape[0]):
        for j in range(a.u.shape[1]):
            Q[i, j] = -0.5*(a.derivative['dudx'][i, j]**2 + a.derivative['dvdy'][i, j]**2) \
            - a.derivative['dudy'][i, j] * a.derivative['dvdx'][i, j]
    return Q

def delta_criterion(a):
    """
    Delta Criterion
    """
    print("Detection method: Delta criterion")
    Q = np.zeros((a.u[0].size, a.u[0].size))
    R = np.zeros((a.u[0].size, a.u[0].size))
    delta = np.zeros((a.u[0].size, a.u[0].size))
    for i in range(a.u[0].size):
        for j in range(a.u[0].size):
            Q[i, j] = -0.5*(a.derivative['dudx'][i, j]**2 + a.derivative['dvdy'][i, j]**2) \
            - a.derivative['dudy'][i, j] * a.derivative['dvdx'][i, j]
            R[i, j] = a.derivative['dudx'][i, j]*a.derivative['dvdy'][i, j] \
            - a.derivative['dvdx'][i, j]*a.derivative['dudy'][i, j]
            delta[i, j] = (Q[i, j]/3)**3 + (R[i, j]/2)**2
    return delta
