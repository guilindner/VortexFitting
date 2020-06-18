#!/usr/bin/env/ python3
"""
Difference schemes for the velocity derivatives
"""

import numpy as np


def second_order_diff(vfield):
    """
    Second order accurate finite difference scheme.

    .. note:: Scheme 1:0:-1

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: class VelocityField
    :returns: derivative
    :rtype: 2D array of float
    """
    print('Difference scheme: Second Order Scheme')

    dx = vfield.x_coordinate_step  # only for homogeneous mesh
    dy = vfield.y_coordinate_step  # only for homogeneous mesh

    vfield.derivative['dudy'], vfield.derivative['dudx'] = np.gradient(vfield.u_velocity_matrix, dy, dx)
    vfield.derivative['dvdy'], vfield.derivative['dvdx'] = np.gradient(vfield.v_velocity_matrix, dy, dx)

    return vfield.derivative


def least_square_diff(vfield):  # may have a problem on the boundary
    """
    Least-square filter difference scheme [RAFFEL1998]_.

    .. note:: Scheme -2:-1:0:1:2

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: class VelocityField
    :returns: derivative
    :rtype: 2D array of float

    .. [RAFFEL1998] Raffel M, Willert C and Kompenhans J 1998
           *Particle image velocimetry A Practical Guide* (Berlin: Springer)
    """
    print('Difference scheme: least-square filter')

    dx = vfield.x_coordinate_step  # only for homogeneous mesh
    dy = vfield.y_coordinate_step  # only for homogeneous mesh

    vfield.derivative['dudx'][:, 2:-2] = (-2 * vfield.u_velocity_matrix[:, 0:-4] -
                                          vfield.u_velocity_matrix[:, 1:-3] +
                                          vfield.u_velocity_matrix[:, 3:-1] +
                                          2 * vfield.u_velocity_matrix[:, 4:]) / (10 * dy)
    vfield.derivative['dudy'][2:-2, :] = (-2 * vfield.u_velocity_matrix[0:-4, :] -
                                          vfield.u_velocity_matrix[1:-3, :] +
                                          vfield.u_velocity_matrix[3:-1, :] +
                                          2 * vfield.u_velocity_matrix[4:, :]) / (10 * dx)
    vfield.derivative['dvdx'][:, 2:-2] = (-2 * vfield.v_velocity_matrix[:, 0:-4] -
                                          vfield.v_velocity_matrix[:, 1:-3] +
                                          vfield.v_velocity_matrix[:, 3:-1] +
                                          2 * vfield.v_velocity_matrix[:, 4:]) / (10 * dy)
    vfield.derivative['dvdy'][2:-2, :] = (-2 * vfield.v_velocity_matrix[0:-4, :] -
                                          vfield.v_velocity_matrix[1:-3, :] +
                                          vfield.v_velocity_matrix[3:-1, :] +
                                          2 * vfield.v_velocity_matrix[4:, :]) / (10 * dx)

    return vfield.derivative


def fourth_order_diff(vfield):
    """
    Fourth order accurate finite difference scheme.

    .. note:: Scheme: 1:-8:0:8:-1

    :params vfield: 2D array of the velocity field, containing u and v
    :type vfield: class VelocityField
    :returns: derivative
    :rtype: 2D array of float
    """
    print('Difference scheme: Fourth Order Scheme')

    dx = vfield.x_coordinate_step  # only for homogeneous mesh
    dy = vfield.y_coordinate_step  # only for homogeneous mesh

    vfield.derivative['dudx'][:, 2:-2] = (vfield.u_velocity_matrix[:, 0:-4] -
                                          8 * vfield.u_velocity_matrix[:, 1:-3] +
                                          8 * vfield.u_velocity_matrix[:, 3:-1] -
                                          vfield.u_velocity_matrix[:, 4:]) / (12 * dy)
    vfield.derivative['dudy'][2:-2, :] = (vfield.u_velocity_matrix[0:-4, :] -
                                          8 * vfield.u_velocity_matrix[1:-3, :] +
                                          8 * vfield.u_velocity_matrix[3:-1, :] -
                                          vfield.u_velocity_matrix[4:, :]) / (12 * dx)
    vfield.derivative['dvdx'][:, 2:-2] = (vfield.v_velocity_matrix[:, 0:-4] -
                                          8 * vfield.v_velocity_matrix[:, 1:-3] +
                                          8 * vfield.v_velocity_matrix[:, 3:-1] -
                                          vfield.v_velocity_matrix[:, 4:]) / (12 * dy)
    vfield.derivative['dvdy'][2:-2, :] = (vfield.v_velocity_matrix[0:-4, :] -
                                          8 * vfield.v_velocity_matrix[1:-3, :] +
                                          8 * vfield.v_velocity_matrix[3:-1, :] -
                                          vfield.v_velocity_matrix[4:, :]) / (12 * dx)

    return vfield.derivative
