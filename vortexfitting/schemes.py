#!/usr/bin/env python3
"""
Difference schemes for the velocity derivatives (updated with SciPy convolutions).

.. tip:: new difference schemes should be added here
"""

import numpy as np
from scipy.ndimage import convolve
from classes import VelocityField


def second_order_diff(vfield: VelocityField) -> np.ndarray:
    """
    Second order accurate finite difference scheme using SciPy convolutions.

    .. note:: Kernel: [-1, 0, 1] / (2dx)

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: class VelocityField

    :returns: derivative
    :rtype: 2D array of float
    """
    print("Difference scheme: Second Order (SciPy convolution)")

    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step

    # Kernels 1D
    kx = np.array([[ -1, 0, 1 ]]) / (2 * dx)   # derivative along x
    ky = np.array([[ -1], [ 0], [ 1 ]]) / (2 * dy)   # derivative along y

    u = vfield.u_velocity_matrix
    v = vfield.v_velocity_matrix

    # Compute derivatives
    vfield.derivative['dudx'] = convolve(u, kx, mode="nearest")
    vfield.derivative['dudy'] = convolve(u, ky, mode="nearest")
    vfield.derivative['dvdx'] = convolve(v, kx, mode="nearest")
    vfield.derivative['dvdy'] = convolve(v, ky, mode="nearest")

    return vfield.derivative



def fourth_order_diff(vfield: VelocityField) -> np.ndarray:
    """
    Fourth order accurate finite difference scheme using SciPy convolutions.

    .. note:: Kernel: [1, -8, 0, 8, -1] / (12 dx)

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: class VelocityField

    :returns: derivative
    :rtype: 2D array of float
    """
    print("Difference scheme: Fourth Order (SciPy convolution)")

    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step

    # Kernels
    kx = np.array([[ 1, -8, 0, 8, -1 ]]) / (12 * dx)
    ky = np.array([[ 1], [-8], [ 0], [ 8], [-1] ]) / (12 * dy)

    u = vfield.u_velocity_matrix
    v = vfield.v_velocity_matrix

    vfield.derivative['dudx'] = convolve(u, kx, mode="nearest")
    vfield.derivative['dudy'] = convolve(u, ky, mode="nearest")
    vfield.derivative['dvdx'] = convolve(v, kx, mode="nearest")
    vfield.derivative['dvdy'] = convolve(v, ky, mode="nearest")

    return vfield.derivative



def least_square_diff(vfield: VelocityField) -> np.ndarray:
    """
    Least-square 5-point derivative scheme via SciPy convolution.

    .. note:: Kernel (Raffel 1998):    [-2, -1, 0, 1, 2] / (10 dx)

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: class VelocityField

    :returns: derivative
    :rtype: 2D array of float
    """
    print("Difference scheme: Least-square (SciPy convolution)")

    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step

    # Raffel LSQ kernel
    kx = np.array([[-2, -1, 0, 1, 2]]) / (10 * dx)
    ky = np.array([[-2], [-1], [ 0], [ 1], [ 2]]) / (10 * dy)

    u = vfield.u_velocity_matrix
    v = vfield.v_velocity_matrix

    vfield.derivative['dudx'] = convolve(u, kx, mode="nearest")
    vfield.derivative['dudy'] = convolve(u, ky, mode="nearest")
    vfield.derivative['dvdx'] = convolve(v, kx, mode="nearest")
    vfield.derivative['dvdy'] = convolve(v, ky, mode="nearest")

    return vfield.derivative