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
    :returns: swirling strength criterion
    :rtype: ndarray
    """
    print('Detection method: 2D swirling strength')

    matrix_a = np.array(
        [[vfield.derivative['dudx'].ravel(), vfield.derivative['dudy'].ravel(),
          vfield.derivative['dudz'].ravel()],
         [vfield.derivative['dvdx'].ravel(), vfield.derivative['dvdy'].ravel(),
          vfield.derivative['dvdz'].ravel()],
         [vfield.derivative['dwdx'].ravel(), vfield.derivative['dwdy'].ravel(),
          -vfield.derivative['dudx'].ravel() - vfield.derivative['dvdy'].ravel()]])

    matrix_a = matrix_a.transpose((2, 1, 0))
    eigenvalues = np.linalg.eigvals(matrix_a)
    swirling = np.max(eigenvalues.imag, axis=1).reshape(vfield.u_velocity_matrix[:, 0].size,
                                                        vfield.u_velocity_matrix[0, :].size)
    print('Max value of swirling: ', np.round(np.max(swirling), 2))
    return swirling


def calc_q_criterion(vfield):
    """
    Q Criterion
    vorticity magnitude and mean strain rate

    :param vfield: contains spatial mesh and velocity components
    :type  vfield: class VelocityField()
    :returns: Q criterion
    :rtype: ndarray
    """
    print('Detection method: Q criterion')
    q_matrix = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    for i in range(vfield.u_velocity_matrix.shape[0]):
        for j in range(vfield.u_velocity_matrix.shape[1]):
            q_matrix[i, j] = -0.5 * (vfield.derivative['dudx'][i, j] ** 2 + vfield.derivative['dvdy'][i, j] ** 2) \
                             - vfield.derivative['dudy'][i, j] * vfield.derivative['dvdx'][i, j]
    return q_matrix


def calc_delta_criterion(vfield):
    """
    Delta Criterion

    :param vfield: contains spatial mesh and velocity components
    :type  vfield: class VelocityField
    :returns: delta criterion
    :rtype: ndarray
    """
    print('Detection method: Delta criterion')
    q_matrix = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    r_matrix = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    delta = np.zeros((vfield.u_velocity_matrix.shape[0], vfield.u_velocity_matrix.shape[1]))
    for i in range(vfield.u_velocity_matrix.shape[0]):
        for j in range(vfield.u_velocity_matrix.shape[1]):
            q_matrix[i, j] = -0.5 * (vfield.derivative['dudx'][i, j] ** 2 + vfield.derivative['dvdy'][i, j] ** 2) \
                             - vfield.derivative['dudy'][i, j] * vfield.derivative['dvdx'][i, j]
            r_matrix[i, j] = vfield.derivative['dudx'][i, j] * vfield.derivative['dvdy'][i, j] \
                             - vfield.derivative['dvdx'][i, j] * vfield.derivative['dudy'][i, j]  # noqa: E261
            delta[i, j] = (q_matrix[i, j] / 3) ** 3 + (r_matrix[i, j] / 2) ** 2
    return delta
