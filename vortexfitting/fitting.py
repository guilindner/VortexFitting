#!/usr/bin/env/ python3
"""
Different functions for the fitting of vortices
"""

import sys
import os
import numpy as np
import scipy.ndimage
import scipy.optimize as opt
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import time
from typing import Any, Optional

sys.path.insert(0, os.path.abspath('./vortexfitting'))
from classes import VelocityField

np.seterr(divide='ignore', invalid='ignore')


def get_fluctuations(velocity_matrix: np.ndarray, mean: np.ndarray, homogeneous_axis: str) -> np.ndarray:
    """
    Used when you have an advective velocity along one axis

    :param velocity_matrix: velocity field
    :type velocity_matrix: ndarray
    :param mean: advective velocity to subtract
    :type mean: np.ndarray
    :param homogeneous_axis: False, 'x', or 'y'. The axis which the mean is subtracted
    :type homogeneous_axis: str

    :returns: input array, minus the advective velocity
    :rtype: ndarray
    """
    if homogeneous_axis is None:
        velocity_matrix = velocity_matrix - mean
    elif homogeneous_axis == 'x':
        velocity_matrix = velocity_matrix - mean[:, None]
    elif homogeneous_axis == 'y':
        velocity_matrix = velocity_matrix - mean[None, :]
    else:
        sys.exit("Invalid homogeneity axis on get_fluctuations function.")
    return velocity_matrix


def normalize(velocity_matrix: np.ndarray, homogeneous_axis: str) -> np.ndarray:
    """
    Normalize with swirling strength

    :param velocity_matrix: velocity field
    :type velocity_matrix: ndarray
    :param homogeneous_axis: False, 'x', or 'y'. The axis which the mean is subtracted
    :type homogeneous_axis: str

    :returns: normalized array
    :rtype: ndarray
    """
    if homogeneous_axis is None:
        velocity_matrix = velocity_matrix / np.sqrt(np.mean(velocity_matrix**2))
    elif homogeneous_axis == 'x':
        velocity_matrix = velocity_matrix / np.sqrt(np.mean(velocity_matrix**2, axis=1))
    elif homogeneous_axis == 'y':
        velocity_matrix = velocity_matrix / np.sqrt(np.mean(velocity_matrix**2, axis=0))
    else:
        sys.exit('Invalid homogeneity axis on normalize function.')
    return velocity_matrix


def window(vfield: VelocityField, x_center_index: int, y_center_index: int, dist: int, theoretical_model: str):
    """
    Defines a window around (x; y) coordinates

    :param theoretical_model: chosen model
    :type theoretical_model: str
    :param vfield: full size velocity field
    :type vfield: VelocityField
    :param x_center_index: box center index (x)
    :type x_center_index: int
    :param y_center_index: box center index (y)
    :type y_center_index: int
    :param dist: size of the vortex (mesh units)
    :param dist: int

    :returns: cropped arrays for x, y, u and v
    :rtype: 2D arrays of floats

    """
    if x_center_index - dist > 0:
        x1 = x_center_index - dist
    else:
        x1 = 0
    if y_center_index - dist > 0:
        y1 = y_center_index - dist
    else:
        y1 = 0
    if x_center_index + dist <= vfield.u_velocity_matrix.shape[1]:
        x2 = x_center_index + dist
    else:
        x2 = vfield.u_velocity_matrix.shape[1]
    if y_center_index + dist <= vfield.v_velocity_matrix.shape[0]:
        y2 = y_center_index + dist
    else:
        y2 = vfield.v_velocity_matrix.shape[0]
    x_index, y_index = np.meshgrid(
        vfield.x_coordinate_matrix[int(x1) : int(x2)], vfield.y_coordinate_matrix[int(y1) : int(y2)], indexing='xy'
    )

    u_data = vfield.u_velocity_matrix[int(y1) : int(y2), int(x1) : int(x2)]
    v_data = vfield.v_velocity_matrix[int(y1) : int(y2), int(x1) : int(x2)]
    if theoretical_model == 'batchelor':
        w_data = vfield.w_velocity_matrix[int(y1) : int(y2), int(x1) : int(x2)]
        return x_index, y_index, u_data, v_data, w_data
    else:
        return x_index, y_index, u_data, v_data


def find_peaks(data: np.ndarray, threshold: float, box_size: int) -> tuple[Any, Any, np.ndarray] | tuple[Any, ...]:
    """
    Find local peaks in an image that are above a specified
    threshold value.

    Peaks are the maxima above the "threshold" within a local region.
    The regions are defined by the "box_size" parameters.
    "box_size" defines the local region around each pixel
    as a square box.

    :param data: The 2D array of the image/data.
    :param threshold: The data value or pixel-wise data values to be used for the
        detection threshold.  A 2D "threshold" must have the same
        shape as "data".
    :param box_size: The size of the local region to search for peaks at every point
    :type data: ndarray
    :type threshold: float
    :type box_size: int

    :returns: An array containing the x and y pixel location of the peaks and their values.
    :rtype: tuple
    """

    if np.all(data == data.flat[0]):
        return tuple()

    data_max = scipy.ndimage.maximum_filter(data, size=box_size, mode='constant', cval=0.0)

    peak_goodmask = data == data_max  # good pixels are True

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[y_peaks, x_peaks]
    peaks = (y_peaks, x_peaks, peak_values)
    return peaks


def direction_rotation(vorticity: np.ndarray, peaks: tuple[Any, Any, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """
    Identify the direction of the vortices rotation using the vorticity.

    :param vorticity: 2D array with the computed vorticity
    :param peaks: list of the detected peaks
    :type vorticity: ndarray
    :type peaks: list

    :returns: vortices_clockwise, vortices_counterclockwise, arrays containing the direction of rotation for each vortex
    :rtype: list
    """

    vortices_clockwise_x, vortices_clockwise_y, vortices_clockwise_cpt = [], [], []
    vortices_counterclockwise_x, vortices_counterclockwise_y, vortices_counterclockwise_cpt = [], [], []
    for i in range(len(peaks[0])):
        if vorticity[peaks[0][i], peaks[1][i]] > 0.0:
            vortices_clockwise_x.append(peaks[0][i])
            vortices_clockwise_y.append(peaks[1][i])
            vortices_clockwise_cpt.append(peaks[2][i])
        else:
            vortices_counterclockwise_x.append(peaks[0][i])
            vortices_counterclockwise_y.append(peaks[1][i])
            vortices_counterclockwise_cpt.append(peaks[2][i])
    vortices_clockwise = (vortices_clockwise_x, vortices_clockwise_y, vortices_clockwise_cpt)
    vortices_counterclockwise = (
        vortices_counterclockwise_x,
        vortices_counterclockwise_y,
        vortices_counterclockwise_cpt,
    )
    vortices_clockwise = np.asarray(vortices_clockwise)
    vortices_counterclockwise = np.asarray(vortices_counterclockwise)
    return vortices_clockwise, vortices_counterclockwise


def correlation_coef(
    u_data: np.ndarray,
    v_data: np.ndarray,
    w_data: Any | None,
    u_model: np.ndarray,
    v_model: np.ndarray,
    w_model: Any | None,
) -> float:
    """Calculates the correlation coefficient between two 2D arrays

    :param u_data: velocity u from the data at the proposed window
    :param v_data: velocity v from the data at the proposed window
    :param w_data: velocity w from the data at the proposed window
    :param u_model: velocity u from the calculated model
    :param v_model: velocity v from the calculated model
    :param w_model: velocity w from the calculated model
    :type u_data: ndarray
    :type v_data: ndarray
    :type w_data: ndarray
    :type u_model: ndarray
    :type v_model: ndarray
    :type w_model: ndarray

    :returns: correlation
    :rtype: float
    """

    u_data = u_data.ravel()
    v_data = v_data.ravel()

    u = u_model.ravel()
    v = v_model.ravel()
    if (w_data is not None) and (w_model is not None):
        w_data = w_data.ravel()
        w = w_model.ravel()

        prod_piv_mod = np.mean(u_data * u + v_data * v + w_data * w)
        prod_piv = np.mean(u * u + v * v + w * w)
        prod_mod = np.mean(u_data * u_data + v_data * v_data + w_data * w_data)
    else:
        prod_piv_mod = np.mean(u_data * u + v_data * v)
        prod_piv = np.mean(u * u + v * v)
        prod_mod = np.mean(u_data * u_data + v_data * v_data)

    correlation = float(prod_piv_mod) / float(max(prod_piv, prod_mod))

    return correlation


def velocity_model(
    core_radius: float,
    gamma: float,
    x_real: float,
    y_real: float,
    u_advection: float,
    v_advection: float,
    x: float,
    y: float,
    uz0: float | None,
    theoretical_model: str,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generates the vortex velocity array, according to the chosen model

    .. tip:: new theoretical should be added here

    :param core_radius: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x_real: relative x position of the vortex center
    :param y_real: relative y position of the vortex center
    :param u_advection: u advective velocity at the center
    :param v_advection: v advective velocity at the center
    :param uz0: vertical velocity
    :param theoretical_model: chosen model
    :param x: x position
    :param y: y position
    :type gamma: float
    :type x_real: float
    :type y_real: float
    :type u_advection: float
    :type v_advection: float
    :type x: float
    :type y: float
    :type uz0: float
    :type theoretical_model: str
    :returns: velx, vely
    :rtype: np.ndarray
    """
    r = np.hypot(x - x_real, y - y_real)

    if theoretical_model == 'lamb-oseen':
        velocity_theta = (gamma / (2 * np.pi * r)) * (1 - np.exp(-(r**2) / core_radius**2))
        velocity_theta = np.nan_to_num(velocity_theta)
        velx = u_advection - velocity_theta * (y - y_real) / r
        vely = v_advection + velocity_theta * (x - x_real) / r
        return np.asarray(np.nan_to_num(velx)), np.asarray(np.nan_to_num(vely))

    elif theoretical_model == 'rankine':
        mask = r < core_radius
        velocity_theta = np.zeros_like(r)

        velocity_theta[mask] = (gamma * r[mask]) / (2 * np.pi * core_radius**2)
        velocity_theta[~mask] = gamma / (2 * np.pi * r[~mask])
        velocity_theta = np.nan_to_num(velocity_theta)
        velx = u_advection - velocity_theta * (y - y_real) / r
        vely = v_advection + velocity_theta * (x - x_real) / r
        return np.nan_to_num(velx), np.nan_to_num(vely)

    elif theoretical_model == 'batchelor':
        velocity_theta = (gamma / (2 * np.pi * r)) * (1 - np.exp(-(r**2) / core_radius**2))
        velocity_theta = np.nan_to_num(velocity_theta)
        velx = u_advection - velocity_theta * (y - y_real) / r
        vely = v_advection + velocity_theta * (x - x_real) / r
        if uz0 is not None:
            velz = uz0 * np.exp(-(r**2) / core_radius**2)
            return np.nan_to_num(velx), np.nan_to_num(vely), np.nan_to_num(velz)
        else:
            return np.nan_to_num(velx), np.nan_to_num(vely)
    else:
        return np.zeros_like(r), np.zeros_like(r)

def process_vortex(i, peaks, vfield, vorticity, rmax, correlation_threshold, theoretical_model, dx, dy):
    """
    Process a single detected vortex (used in parallel execution).

    :param i: index of the vortex to process
    :type i: int
    :param peaks: tuple of arrays containing coordinates of detected peaks (y_indices, x_indices)
    :type peaks: tuple
    :param vfield: velocity field object containing u, v and metadata
    :type vfield: VelocityField
    :param vorticity: 2D array of vorticity values
    :type vorticity: np.ndarray
    :param rmax: maximum core radius; if equal to 0.0, a default radius based on dx and dy is used
    :type rmax: float
    :param correlation_threshold: minimum correlation coefficient required to accept the vortex
    :type correlation_threshold: float
    :param theoretical_model: vortex model to fit ('lamb-oseen', 'rankine', or 'batchelor')
    :type theoretical_model: str
    :param dx: grid spacing in the x direction
    :type dx: float
    :param dy: grid spacing in the y direction
    :type dy: float

    :returns: a list of fitted vortex parameters if the vortex is valid and passes the
              correlation threshold, or ``None`` otherwise.
              Returned format:

              For Lamb–Oseen / Rankine models::
                  [core_radius, gamma, x_center, y_center,
                   u_adv, v_adv, window_radius, correlation,
                   u_theta, 0.0]

              For Batchelor model::
                  [core_radius, gamma, x_center, y_center,
                   u_adv, v_adv, window_radius, correlation,
                   u_theta, uz0]

    :rtype: list or None
    """

    x_center_index = peaks[1][i]
    y_center_index = peaks[0][i]
    # print(i, 'Processing detected swirling at (x, y)', x_center_index, y_center_index)

    core_radius = 2 * np.hypot(dx, dy) if rmax == 0.0 else rmax
    gamma = vorticity[y_center_index, x_center_index] * np.pi * core_radius**2
    velz0 = np.max(np.hypot(vfield.u_velocity_matrix, vfield.v_velocity_matrix))

    # Fit vortex parameters
    vortices_parameters = full_fit(core_radius, gamma, vfield, x_center_index, y_center_index, velz0, theoretical_model)

    if vortices_parameters[6] < 2:
        return None  # Skip invalid vortices

    # Extract velocity data
    if theoretical_model == 'batchelor':
        x_index, y_index, u_data, v_data, w_data = window(
            vfield,
            round(vortices_parameters[2] / dx),
            round(vortices_parameters[3] / dy),
            int(vortices_parameters[6]),
            theoretical_model,
        )
        u_model, v_model, w_model = velocity_model(
            vortices_parameters[0],
            vortices_parameters[1],
            vortices_parameters[2],
            vortices_parameters[3],
            vortices_parameters[4],
            vortices_parameters[5],
            x_index,
            y_index,
            np.mean(w_data),
            theoretical_model,
        )
    else:
        x_index, y_index, u_data, v_data = window(
            vfield,
            round(vortices_parameters[2] / dx),
            round(vortices_parameters[3] / dy),
            int(vortices_parameters[6]),
            theoretical_model,
        )
        u_model, v_model = velocity_model(
            vortices_parameters[0],
            vortices_parameters[1],
            vortices_parameters[2],
            vortices_parameters[3],
            vortices_parameters[4],
            vortices_parameters[5],
            x_index,
            y_index,
            None,
            theoretical_model,
        )

    # Compute correlation coefficient
    correlation_value = correlation_coef(
        u_data - vortices_parameters[4],
        v_data - vortices_parameters[5],
        None,
        u_model - vortices_parameters[4],
        v_model - vortices_parameters[5],
        None,
    )

    if correlation_value > correlation_threshold:
        print(f'Accepted! Correlation = {correlation_value:.2f} (vortex #{i})')
        u_theta = (vortices_parameters[1] / (2 * np.pi * vortices_parameters[0])) * (1 - np.exp(-1))

        if theoretical_model == 'batchelor':
            return [
                vortices_parameters[0],
                vortices_parameters[1],
                vortices_parameters[2],
                vortices_parameters[3],
                vortices_parameters[4],
                vortices_parameters[5],
                vortices_parameters[6],
                correlation_value,
                u_theta,
                vortices_parameters[7],
            ]
        else:
            return [
                vortices_parameters[0],
                vortices_parameters[1],
                vortices_parameters[2],
                vortices_parameters[3],
                vortices_parameters[4],
                vortices_parameters[5],
                vortices_parameters[6],
                correlation_value,
                u_theta,
                0.0,
            ]
    return None  # Skip non-matching vortices


def get_vortices(
    vfield: VelocityField,
    peaks: tuple[Any, Any, np.ndarray],
    vorticity: np.ndarray,
    rmax: float,
    correlation_threshold: float,
    theoretical_model: str,
    num_cores: Optional[int] = 1,
) -> list:
    """
    General routine to check if the detected vortex is a real vortex.
    Parallelized using joblib.

    :returns: list of detected vortices
    """
    start_time = time.time()  # Start timing

    dx, dy = vfield.x_coordinate_step, vfield.y_coordinate_step
    print(f"Using {num_cores} cores") 
    results = Parallel(n_jobs=num_cores,
                       backend="threading")(
        delayed(process_vortex)(i, peaks, vfield, vorticity, rmax, correlation_threshold, theoretical_model, dx, dy)
        for i in range(len(peaks[0]))
    )

    # Filter out None values
    vortices = [v for v in results if v is not None]

    end_time = time.time()  # End timing
    print(f"Time taken for get_vortices: {end_time - start_time:.3f} seconds")

    return vortices


def full_fit(
    core_radius: float,
    gamma: float,
    vfield: VelocityField,
    x_center_index: int,
    y_center_index: int,
    velz0: float,
    theoretical_model: str,
) -> (
    tuple[float, float, float, float, float, float, float, float]
    | tuple[float, float, float, float, float, float, float]
):
    """Full fitting procedure

    :param velz0: vertical velocity at the vortex center (Batchelor's model)
    :type velz0: float
    :param core_radius: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param vfield: data from the input file
    :param x_center_index: x index of the vortex center
    :param y_center_index: y index of the vortex center
    :type core_radius: float
    :type gamma: float
    :type vfield: class
    :type x_center_index: int
    :type y_center_index: int
    :param theoretical_model: chosen model
    :type theoretical_model: str
    :returns: fitted[i], dist
    :rtype: list
    """

    fitted = [0.0]*7
    fitted[0] = core_radius
    fitted[1] = gamma
    fitted[2] = vfield.x_coordinate_matrix[x_center_index]
    fitted[3] = vfield.y_coordinate_matrix[y_center_index]
    if theoretical_model == 'batchelor':
        fitted[6] = velz0
    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step
    dist = 0
    # correlation_value = 0.0
    for i in range(10):

        x_center_index = int(round(fitted[2] / dx))
        y_center_index = int(round(fitted[3] / dy))
        if x_center_index >= vfield.u_velocity_matrix.shape[1]:
            x_center_index = vfield.u_velocity_matrix.shape[1] - 1
        if x_center_index <= 2:
            x_center_index = 3
        if y_center_index >= vfield.v_velocity_matrix.shape[0]:
            y_center_index = vfield.v_velocity_matrix.shape[0] - 1
        r1 = fitted[0]
        x1 = fitted[2]
        y1 = fitted[3]
        if theoretical_model == 'batchelor':
            velz0 = fitted[6]
        dist = int(round(fitted[0] / np.hypot(dx, dy), 0)) + 1
        if fitted[0] < 2 * np.hypot(dx, dy):
            break
        fitted[4] = vfield.u_velocity_matrix[y_center_index, x_center_index]  # u_advection
        fitted[5] = vfield.v_velocity_matrix[y_center_index, x_center_index]  # v_advection

        if theoretical_model == 'batchelor':
            x_index, y_index, u_data, v_data, w_data = window(
                vfield, x_center_index, y_center_index, dist, theoretical_model
            )
            fitted = fit(
                fitted[0],
                fitted[1],
                x_index,
                y_index,
                fitted[2],
                fitted[3],
                u_data,
                v_data,
                w_data,
                fitted[6],
                fitted[4],
                fitted[5],
                i,
                theoretical_model,
            )
        else:
            x_index, y_index, u_data, v_data = window(vfield, x_center_index, y_center_index, dist, theoretical_model)
            fitted = fit(
                fitted[0],
                fitted[1],
                x_index,
                y_index,
                fitted[2],
                fitted[3],
                u_data,
                v_data,
                None,
                None,
                fitted[4],
                fitted[5],
                i,
                theoretical_model,
            )

        if i > 0:
            # break if radius variation is less than 10% and accepts
            if abs(fitted[0] / r1 - 1) < 0.1:
                if (abs((fitted[2] / x1 - 1)) < 0.1) or (abs((fitted[3] / y1 - 1)) < 0.1):
                    break
            # break if x or y position is out of the window and discards
            if (abs((fitted[2] - x1)) > dist * dx) or (abs((fitted[3] - y1)) > dist * dy):
                dist = 0
                break
    # return core_radius, gamma, xcenter, ycenter, u_advection, v_advection, dist
    if theoretical_model == 'batchelor':
        return fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], dist, fitted[6]

    else:
        return fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], dist


def fit(
    core_radius: float,
    gamma: float,
    x: np.ndarray,
    y: np.ndarray,
    x_real: float,
    y_real: float,
    u_data: np.ndarray,
    v_data: np.ndarray,
    w_data: Any | None,
    velz0: float | None,
    u_advection: float,
    v_advection: float,
    i: int,
    theoretical_model: str,
) -> np.ndarray:
    """
    Fitting of the vortex with the desired theoretical model

    .. tip:: new theoretical should be added here

    :param w_data: vertical velocity
    :type w_data:  np.ndarray
    :param velz0: vertical velocity at the vortex center (Batchelor's model)
    :type velz0: float
    :param core_radius: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x: x position
    :param y: y position
    :param x_real: x position of the vortex center
    :param y_real: y position of the vortex center
    :param u_data: velocity u from the data at the proposed window
    :param v_data: velocity v from the data at the proposed window
    :param u_advection: uniform advection velocity u
    :param v_advection: uniform advection velocity u
    :param i: current iteration for fitting
    :type core_radius: float
    :type gamma: float
    :type x: ndarray
    :type y: ndarray
    :type x_real: float
    :type y_real: float
    :type u_data: ndarray
    :type v_data: ndarray
    :type u_advection: float
    :type v_advection: float
    :type i: iterator
    :param theoretical_model: chosen model
    :type theoretical_model: str

    :returns: fitted parameters (core_radius, gamma,xcenter,ycenter, u_advection, v_advection...)
    :rtype: list
    """
    # Method for opt.least_squares fitting. Can be
    # 'trf': Trust Region Reflective algorithm
    # 'dogbox': dogleg algorithm
    # 'lm': Levenberg-Marquardt algorithm
    method = 'trf'
    epsilon_r = 1e-10
    epsilon_model = 0.001
    applied_model = ''
    x = x.ravel()
    y = y.ravel()
    u_data = u_data.ravel()
    v_data = v_data.ravel()
    
    if theoretical_model == 'batchelor':
        w_data = w_data.ravel()
    else:
        w_data = np.zeros_like(u_data)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    def rankine_model(fitted: list) -> np.ndarray:
        """
        Rankine velocity model used for the nonlinear fitting

        :param fitted: parameters of a vortex (core_radius, gamma,xcenter,ycenter, u_advection, v_advection)
        :type fitted: list
        :returns: velocity field, following a Rankine model
        :rtype: ndarray
        """

        core_radius_model = fitted[0]
        gamma_model = fitted[1]
        xcenter_model = fitted[2]
        ycenter_model = fitted[3]
        u_advection_model = fitted[4]
        v_advection_model = fitted[5]

        r = np.hypot(x - xcenter_model, y - ycenter_model) + epsilon_r

        mask = r < core_radius
        u_theta_model = np.zeros_like(r)
        u_theta_model[mask] = (gamma_model * r[mask]) / (2 * np.pi * core_radius_model**2)
        u_theta_model[~mask] = gamma_model / (2 * np.pi * r[~mask])
        u_theta_model = np.nan_to_num(u_theta_model)
        u_model = u_advection_model - u_theta_model * (y - ycenter_model) / r - u_data
        v_model = v_advection_model + u_theta_model * (x - xcenter_model) / r - v_data
        u_model = np.nan_to_num(u_model)
        v_model = np.nan_to_num(v_model)
        # w_model = np.zeros_like(u_model)
        vfield_model = np.append(u_model, v_model)
        return vfield_model

    def lamb_oseen_model(fitted: list) -> np.ndarray:
        """
        Lamb-Oseen velocity model used for the nonlinear fitting

        :param fitted: parameters of a vortex (core_radius, gamma,xcenter,ycenter, u_advection, v_advection)
        :type fitted: list
        :returns: velocity field, following a Lamb-Oseen model
        :rtype: ndarray
        """

        core_radius_model = fitted[0]
        gamma_model = fitted[1]
        xcenter_model = fitted[2]
        ycenter_model = fitted[3]
        u_advection_model = fitted[4]
        v_advection_model = fitted[5]
        r = np.hypot(x - xcenter_model, y - ycenter_model) + epsilon_r
        u_theta_model = gamma_model / (2 * np.pi * r) * (1 - np.exp(-(r**2) / core_radius_model**2))
        u_theta_model = np.nan_to_num(u_theta_model)
        u_model = u_advection_model - u_theta_model * (y - ycenter_model) / r - u_data
        v_model = v_advection_model + u_theta_model * (x - xcenter_model) / r - v_data
        u_model = np.nan_to_num(u_model)
        v_model = np.nan_to_num(v_model)
        # w_model = np.zeros_like(u_model)
        vfield_model = np.append(u_model, v_model)
        return vfield_model

    def batchelor_model(fitted: list) -> np.ndarray:
        """
        Batchelor velocity model used for the nonlinear fitting

        :param fitted: parameters of a vortex (core_radius, gamma,xcenter,ycenter, u_advection, v_advection)
        :type fitted: list

        :returns: velocity field, following a Batchelor model
        :rtype: ndarray
        """

        core_radius_model = fitted[0]
        gamma_model = fitted[1]
        xcenter_model = fitted[2]
        ycenter_model = fitted[3]
        u_advection_model = fitted[4]
        v_advection_model = fitted[5]
        velz0_model = fitted[6]
        r = np.hypot(x - xcenter_model, y - ycenter_model) + epsilon_r
        u_theta_model = gamma_model / (2 * np.pi * r) * (1 - np.exp(-(r**2) / core_radius_model**2))
        u_theta_model = np.nan_to_num(u_theta_model)
        u_model = u_advection_model - u_theta_model * (y - ycenter_model) / r - u_data
        v_model = v_advection_model + u_theta_model * (x - xcenter_model) / r - v_data
        w_model = velz0_model * np.exp(-(r**2) / core_radius_model**2) - w_data
        u_model = np.nan_to_num(u_model)
        v_model = np.nan_to_num(v_model)
        w_model = np.nan_to_num(w_model)
        vfield_model = np.concatenate((u_model, v_model, w_model), axis=0)
        return vfield_model

    if i > 0:
        m = 1.0
    else:
        m = 4.0

    match theoretical_model:
        case 'rankine':
            applied_model = rankine_model
        case 'batchelor':
            applied_model = batchelor_model
        case _:
            applied_model = lamb_oseen_model

    if theoretical_model == 'batchelor':
        bnds = (
            [
                0,
                gamma - abs(gamma) * m / 2 - epsilon_model,
                x_real - m * dx - epsilon_model,
                y_real - m * dy - epsilon_model,
                u_advection - abs(u_advection) - epsilon_model,
                v_advection - abs(v_advection) - epsilon_model,
                velz0 - abs(velz0) * m / 2 - epsilon_model,
            ],
            [
                core_radius + core_radius * m,
                gamma + abs(gamma) * m / 2 + epsilon_model,
                x_real + m * dx + epsilon_model,
                y_real + m * dy + epsilon_model,
                u_advection + abs(u_advection) + epsilon_model,
                v_advection + abs(v_advection) + epsilon_model,
                velz0 + abs(velz0) * m / 2 + epsilon_model,
            ],
        )
        if method == 'trf':
            sol = opt.least_squares(
                applied_model,
                [core_radius, gamma, x_real, y_real, u_advection, v_advection, velz0],
                method='trf',
                bounds=bnds,
            )
        elif method == 'dogbox':
            sol = opt.least_squares(
                applied_model,
                [core_radius, gamma, x_real, y_real, u_advection, v_advection, velz0],
                method='dogbox',
                bounds=bnds,
            )
        else:
            sol = opt.least_squares(
                applied_model,
                [core_radius, gamma, x_real, y_real, u_advection, v_advection, velz0],
                method='lm',
                xtol=10 * np.hypot(dx, dy),
            )
    else:
        bnds = (
            [
                0,
                gamma - abs(gamma) * m / 2 - epsilon_model,
                x_real - m * dx - epsilon_model,
                y_real - m * dy - epsilon_model,
                u_advection - abs(u_advection) - epsilon_model,
                v_advection - abs(v_advection) - epsilon_model,
            ],
            [
                core_radius + core_radius * m,
                gamma + abs(gamma) * m / 2 + epsilon_model,
                x_real + m * dx + epsilon_model,
                y_real + m * dy + epsilon_model,
                u_advection + abs(u_advection) + epsilon_model,
                v_advection + abs(v_advection) + epsilon_model,
            ],
        )
        if method == 'trf':
            sol = opt.least_squares(
                applied_model, [core_radius, gamma, x_real, y_real, u_advection, v_advection], method='trf', bounds=bnds
            )
        elif method == 'dogbox':
            sol = opt.least_squares(
                applied_model,
                [core_radius, gamma, x_real, y_real, u_advection, v_advection],
                method='dogbox',
                bounds=bnds,
            )
        else:
            sol = opt.least_squares(
                applied_model,
                [core_radius, gamma, x_real, y_real, u_advection, v_advection],
                method='lm',
                xtol=10 * np.hypot(dx, dy),
            )
            # print(sol.x)
    return sol.x


def plot_fields(vfield: VelocityField, detection_field: np.ndarray) -> None:
    """
    Plot fields: display the (u,v,w) fields and the vorticity field.

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField
    :param detection_field: detection field (vorticity ...)
    :type detection_field: ndarray

    :returns: popup
    :rtype: image
    """
    plt.ion()
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)  # , sharex='col', sharey='row')
    im1 = ax1.imshow(vfield.u_velocity_matrix, cmap='seismic', origin="lower")
    ax1.set_title('Velocity u (velocity_s)')
    fig.colorbar(im1, ax=ax1)
    im2 = ax2.imshow(vfield.v_velocity_matrix, cmap='seismic', origin="lower")
    ax2.set_title('Velocity v (velocity_n)')
    fig.colorbar(im2, ax=ax2)

    try:
        im3 = ax3.imshow(vfield.w_velocity_matrix, cmap='seismic', origin="lower")
        ax3.set_title('Velocity w (velocity_z)')
        fig.colorbar(im3, ax=ax3)
    except AttributeError:
        print('No w velocity')
    im4 = ax4.imshow(detection_field, origin="lower", cmap='seismic')
    ax4.set_title('Vorticity')
    fig.colorbar(im4, ax=ax4)
    plt.tight_layout()

    plt.show()


def plot_detect(
    vortices_counterclockwise: np.ndarray, vortices_clockwise: np.ndarray, detection_field: np.ndarray, flip_axis: bool
) -> None:
    """
    Plot detect: display the location and rotation of the vortices

    :param vortices_counterclockwise: vortices spinning counterclockwise
    :type vortices_counterclockwise: array of vortices
    :param vortices_clockwise: vortices spinning clockwise
    :type vortices_clockwise: array of vortices
    :param detection_field: detection field (vorticity ...)
    :type detection_field: ndarray
    :param flip_axis: for flipping x/y axis
    :type flip_axis: bool

    :returns: popup
    :rtype: image
    """
    # plt.ion()
    plt.subplot()
    if flip_axis:
        detection_field = detection_field.T  # transpose the detection field
        plt.scatter(
            vortices_counterclockwise[0],
            vortices_counterclockwise[1],
            edgecolor='green',
            facecolor='green',
            label='counterclockwise',
        )
        plt.scatter(
            vortices_clockwise[0], vortices_clockwise[1], edgecolor='yellow', facecolor='yellow', label='clockwise'
        )
    else:
        plt.scatter(
            vortices_counterclockwise[1],
            vortices_counterclockwise[0],
            edgecolor='green',
            facecolor='green',
            label='counterclockwise',
        )
        plt.scatter(
            vortices_clockwise[1], vortices_clockwise[0], edgecolor='yellow', facecolor='yellow', label='clockwise'
        )

    plt.title('Detected possible vortices')
    # plt.contourf(field, cmap="Greys_r")

    plt.imshow(detection_field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    # plt.legend()
    # plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()

    plt.show(block=True)


def plot_quiver(
    x_index: np.ndarray, y_index: np.ndarray, u_data: np.ndarray, v_data: np.ndarray, detection_field: np.ndarray
) -> None:
    """
    Plot quiver: display a specific (x,y) location with vector fields.

    :param x_index: contains spatial mesh (x direction)
    :type x_index: ndarray
    :param y_index: contains spatial mesh (y direction)
    :type y_index: ndarray
    :param u_data: contains velocity data (u component)
    :type u_data: ndarray
    :param v_data: contains velocity data (v component)
    :type v_data: ndarray
    :param detection_field: detection field (vorticity ...)
    :type detection_field: ndarray

    :returns: popup
    :rtype: image
    """
    plt.ion()
    plt.figure()
    plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(detection_field, extent=[x_index[0][0], x_index[0][-1], y_index[0][0], y_index[-1][0]])
    s = 1  # sampling factor, can be modified
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s])
    plt.show()


def plot_fit(
    x_index: np.ndarray,
    y_index: np.ndarray,
    u_data: np.ndarray,
    v_data: np.ndarray,
    u_model: np.ndarray,
    v_model: np.ndarray,
    xc: float,
    yc: float,
    core_radius: float,
    gamma: float,
    u_advection: np.ndarray,
    v_advection: np.ndarray,
    correlation_value: float,
    cpt_vortex: int,
    subtract_advection_field: bool,
    output_dir: str,
    time_step: int,
    output_format: str,
) -> None:
    """
    Plot the data velocity field and the model velocity field with a quiver plot

    :param x_index: contains spatial mesh (x direction)
    :type x_index: ndarray
    :param y_index: contains spatial mesh (y direction)
    :type y_index: ndarray
    :param u_data: contains velocity data (u component)
    :type u_data: ndarray
    :param v_data: contains velocity data (v component)
    :type v_data: ndarray
    :param u_model: contains velocity data (u component)
    :type u_model: ndarray
    :param v_model: contains velocity data (v component)
    :type v_model: ndarray
    :param xc: x coordinate of the vortex center
    :type xc: float
    :param yc: y coordinate of the vortex center
    :type yc: float
    :param core_radius: dimension of the vortex core radius
    :type core_radius: float
    :param gamma: circulation of the vortex
    :type gamma: float
    :param u_advection: contains velocity data (u component)
    :type u_advection: ndarray
    :param v_advection: contains velocity data (v component)
    :type v_advection: ndarray
    :param cpt_vortex: current n° of the vortex
    :type cpt_vortex: int
    :param subtract_advection_field: if True, the advection field (u_advection, v_advection) is subtracted
    :type subtract_advection_field:  bool
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param correlation_value: correlation between the vortex and a chosen theoretical model
    :type correlation_value: float
    :param time_step: current time_step
    :type time_step: int
    :param output_format: format for output files (pdf, png ...)
    :type output_format: str

    :returns: image file
    :rtype: image
    """
    plt.figure()
    s = 1  # sampling factor
    if x_index.size > 400:
        s = 2
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s], color='r', label='data')
    plt.quiver(
        x_index[::s, ::s], y_index[::s, ::s], u_model[::s, ::s], v_model[::s, ::s], color='b', label='model', alpha=0.5
    )
    circle1 = plt.Circle((xc, yc), core_radius, color='k', alpha=0.05)
    plt.gca().add_artist(circle1)
    plt.gca().scatter([xc], [yc], marker='+', color='k', s=100)
    plt.legend()
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.title(
        r'r=%s $\Gamma$=%s u=%s v=%s C=%s'
        % (
            round(core_radius, 2),
            round(gamma, 2),
            np.round(u_advection, 2),
            np.round(v_advection, 2),
            round(correlation_value, 2),
        )
    )
    if not subtract_advection_field:
        plt.savefig(
            output_dir + '/vortex%i_%i_%s.%s' % (time_step, cpt_vortex, 'initial_vfield', output_format),
            format=output_format,
        )
    else:
        plt.savefig(
            output_dir + '/vortex%i_%i_%s.%s' % (time_step, cpt_vortex, 'advection_field_subtracted', output_format),
            format=output_format,
        )
    plt.close('all')


def plot_accepted(
    vfield: VelocityField,
    vortices_list: list,
    detection_field: np.ndarray,
    output_dir: str,
    time_step: int,
    output_format: str,
) -> None:
    """
    Plot accepted: display the accepted vortices, with respect to the different criteria
    (correlation threshold, box size ...)

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param detection_field: detection field (vorticity ...)
    :type detection_field: ndarray
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :param output_format: format for output files (pdf, png ...)
    :type output_format: str

    :returns: popup
    :rtype: image
    """
    plt.figure(1)

    plt.contourf(vfield.x_coordinate_matrix, vfield.y_coordinate_matrix, detection_field, origin='lower', cmap="bone")
    plt.xlabel('x')
    plt.ylabel('y')
    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step
    plt.figure(2)
    plt.imshow(detection_field, origin='lower', cmap="bone")
    plt.xlabel('x (mesh units)')
    plt.ylabel('y (mesh units)')
    plt.colorbar()
    for i, line in enumerate(vortices_list):
        if vortices_list[i][1] > 0:
            # Gamma > 0 (vorticity > 0): counterclockwise rotation
            plt.figure(1)
            circle1 = plt.Circle((line[2], line[3]), line[0], edgecolor='green', facecolor='none', gid='vortex%i' % i)
            plt.gca().add_artist(circle1)
            plt.figure(2)
            circle1 = plt.Circle(
                (line[2] / dx, line[3] / dy),
                line[0] / np.hypot(dx, dy),
                edgecolor='green',
                facecolor='none',
                gid='vortex%i' % i,
            )
            plt.gca().add_artist(circle1)
        else:
            # Gamma < 0 (vorticity < 0): clockwise rotation
            plt.figure(1)
            circle1 = plt.Circle((line[2], line[3]), line[0], edgecolor='yellow', facecolor='none', gid='vortex%i' % i)
            plt.gca().add_artist(circle1)
            plt.figure(2)
            circle1 = plt.Circle(
                (line[2] / dx, line[3] / dy),
                line[0] / np.hypot(dx, dy),
                edgecolor='yellow',
                facecolor='none',
                gid='vortex%i' % i,
            )
            plt.gca().add_artist(circle1)

    # Comparing data
    # file_in = open('../data/comparison_data.dat', 'r')
    # for line in file_in:
    #     x_comp = int(float(line.split()[1]))
    #     y_comp = int(float(line.split()[2]))
    #     gamma_comp = float(line.split()[3])
    #     r_comp = float(line.split()[4])
    #     if gamma_comp > 0:
    #         orient = 'green'
    #     else:
    #         orient = 'yellow'
    #     circle2 = plt.Circle((x_comp, y_comp), r_comp, edgecolor=orient, facecolor='none')
    #     plt.gca().add_artist(circle2)

    # plt.legend()
    plt.figure(1)
    plt.tight_layout()
    plt.savefig(output_dir + '/accepted_{}.svg'.format(time_step), format='svg')
    plt.savefig(output_dir + '/accepted_{}.{}'.format(time_step, output_format), format=output_format)
    # create_links(output_dir + '/accepted_{}.svg'.format(time_step), vortices_list, output_dir, time_step,
    #             output_format)
    plt.figure(2)
    plt.savefig(
        output_dir + '/meshed_{:01d}.{}'.format(time_step, output_format), format=output_format
    )  # use it to verify your coordinate system if needed !
    # plt.savefig(output_dir+'/tk_{:01d}.png'.format(time_step), format='png', transparent=True)

    # plt.show()


def plot_vortex(
    vfield: VelocityField,
    vortices_list: list,
    output_dir: str,
    time_step: int,
    output_format: str,
    theoretical_model: str,
) -> None:
    """
    Plot vortex: plot a vortex and its corresponding vortex model

    :param theoretical_model: chosen vortex model for the fitting
    :type theoretical_model: str
    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :param output_format: format for output files (pdf, png ...)
    :type output_format: str

    :returns: file
    :rtype: image
    """
    for cpt_vortex, line in enumerate(vortices_list):
        print(
            'r: {:.3f}'.format(line[0]),
            'gamma: {:.2f}'.format(line[1]),
            'xc: {:.2f}'.format(line[2]),
            'yc: {:.2f}'.format(line[3]),
            'correlation: {:.2f}'.format(line[7]),
            'utheta: {:.2f}'.format(line[8]),
        )
        dx = vfield.x_coordinate_step
        dy = vfield.y_coordinate_step

        if theoretical_model == 'batchelor':
            x_index, y_index, u_data, v_data, w_data = window(
                vfield, round(line[2] / dx, 0), round(line[3] / dy, 0), line[6], theoretical_model
            )
            u_model, v_model, w_model = velocity_model(
                line[0],
                line[1],
                line[2],
                line[3],
                line[4],
                line[5],
                x_index,
                y_index,
                np.mean(w_data),
                theoretical_model,
            )
            print('\tuz0: {:.2f}'.format(line[9]))
        else:
            x_index, y_index, u_data, v_data = window(
                vfield, round(line[2] / dx, 0), round(line[3] / dy, 0), line[6], theoretical_model
            )
            u_model, v_model = velocity_model(
                line[0], line[1], line[2], line[3], line[4], line[5], x_index, y_index, None, theoretical_model
            )
        correlation_value = correlation_coef(u_data, v_data, None, u_model, v_model, None)

        plot_fit(
            x_index,
            y_index,
            u_data,
            v_data,
            u_model,
            v_model,
            line[2],
            line[3],
            line[0],
            line[1],
            line[4],
            line[5],
            correlation_value,
            cpt_vortex,
            False,
            output_dir,
            time_step,
            output_format,
        )
        correlation_value = correlation_coef(
            u_data - line[4], v_data - line[5], None, u_model - line[4], v_model - line[5], None
        )
        plot_fit(
            x_index,
            y_index,
            u_data - line[4],
            v_data - line[5],
            u_model - line[4],
            v_model - line[5],
            line[2],
            line[3],
            line[0],
            line[1],
            line[4],
            line[5],
            correlation_value,
            cpt_vortex,
            True,
            output_dir,
            time_step,
            output_format,
        )


def create_links(path: str, vortices_list: list, output_dir: str, time_step: int, output_format: str) -> None:
    """
    create links: add some links between the accepted.svg file and the detected vortices

    :param path: path of the accepted.svg file
    :type path: str
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :param output_format: format for output files (pdf, png ...)
    :type output_format: str-9888

    :returns: file
    :rtype: image
    """

    file_in = open(path, "r")
    file_out = open(output_dir + "/linked_{:01d}.svg".format(time_step), "w")
    i = 0
    vortex_found = False
    for line in file_in:
        if "</g>" in line:
            if vortex_found:
                file_out.write(line)
                file_out.write('   </a>\n')
                vortex_found = False
            else:
                file_out.write(line)
        elif "vortex" in line:
            file_out.write('   <a href="vortex%i_%i_advection_field_subtracted.%s">\n' % (time_step, i, output_format))
            file_out.write(line)
            file_out.write(
                '   <title>Vortex %i: r = %s gamma = %s</title>\n'
                % (i, round(vortices_list[i][3], 1), round(vortices_list[i][2], 1))
            )
            i += 1
            vortex_found = True
        else:
            file_out.write(line)
