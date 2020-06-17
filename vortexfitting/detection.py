#!/usr/bin/env/ python3
"""
Different methods for vortex detection
"""

import numpy as np


def calc_swirling(vfield):
    """
    2D Swirling strength

    :param vfield: contains spatial mesh and velocity components
    :type  vfield: class VelocityField()
    :returns: swirling
    :rtype: 2D array of float
    """
    print('Detection method: 2D swirling strength')
    A = np.zeros((vfield.u_velocity_matrix.size, 3, 3))
    A = np.array(
        [[vfield.derivative['dudx'].ravel(), vfield.derivative['dudy'].ravel(),
          vfield.derivative['dudz'].ravel()],
         [vfield.derivative['dvdx'].ravel(), vfield.derivative['dvdy'].ravel(),
          vfield.derivative['dvdz'].ravel()],
         [vfield.derivative['dwdx'].ravel(), vfield.derivative['dwdy'].ravel(),
          -vfield.derivative['dudx'].ravel()-vfield.derivative['dvdy'].ravel()]])

    A = A.transpose(2, 1, 0)
    eigenvalues = np.linalg.eigvals(A)
    swirling = np.max(eigenvalues.imag, axis=1).reshape(vfield.u_velocity_matrix[:, 0].size, vfield.u_velocity_matrix[0, :].size)
    print('Max value of swirling: ', np.round(np.max(swirling),2))
    return swirling

def q_criterion(vfield):
    """
    Q Criterion
    vorticity magnitude and mean strain rate

    :param vfield: contains spatial mesh and velocity components
    :type  vfield: class VelocityField()
    :returns: Q
    :rtype: 2D array of float
    """
    print('Detection method: Q criterion')
    Q = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    #print(vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1])
    for i in range(vfield.u_velocity_matrix.shape[0]):
        for j in range(vfield.u_velocity_matrix.shape[1]):
            Q[i, j] = -0.5*(vfield.derivative['dudx'][i, j]**2 + vfield.derivative['dvdy'][i, j]**2) \
            - vfield.derivative['dudy'][i, j] * vfield.derivative['dvdx'][i, j]
    return Q

def delta_criterion(vfield):
    """
    Delta Criterion

    :param vfield: contains spatial mesh and velocity components
    :type  vfield: class VelocityField()
    :returns: delta
    :rtype: 2D array of float
    """
    print('Detection method: Delta criterion')
    Q = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    R = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    delta = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    for i in range(vfield.u_velocity_matrix.shape[0]):
        for j in range(vfield.u_velocity_matrix.shape[1]):
            Q[i, j] = -0.5*(vfield.derivative['dudx'][i, j]**2 + vfield.derivative['dvdy'][i, j]**2) \
            - vfield.derivative['dudy'][i, j] * vfield.derivative['dvdx'][i, j]
            R[i, j] = vfield.derivative['dudx'][i, j]*vfield.derivative['dvdy'][i, j] \
            - vfield.derivative['dvdx'][i, j]*vfield.derivative['dudy'][i, j]
            delta[i, j] = (Q[i, j]/3)**3 + (R[i, j]/2)**2
    return delta
