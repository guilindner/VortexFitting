#!/usr/bin/env/ python3
"""
Different functions for the fitting of vortices
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
np.seterr(divide='ignore', invalid='ignore')
import scipy.optimize as opt


def get_fluc(x, mean, hom_axis):
    """
    Used when you have a advective velocity along one axis

    :param x: velocity field
    :type x: 2D array of float
    :param mean: advective velocity to subtract
    :type mean: float
    :param hom_axis: False, 'x', or 'y'. The axis which the mean is subtracted
    :type hom_axis: str

    :returns: input array, minus the advective velocity
    :rtype: 2D arrays of float
    """
    if hom_axis is None:
        x = x - mean
    elif hom_axis == 'x':
        x = x - mean[:, None]
    elif hom_axis == 'y':
        x = x - mean[None, :]
    else:
        sys.exit("Invalid homogenity axis.")
    return x


def normalize(x, hom_axis):
    """
    Normalize with swirling strength

    :param x: velocity field
    :type x: 2D array of float
    :param hom_axis: False, 'x', or 'y'. The axis which the mean is subtracted
    :type hom_axis: str

    :returns: normalized array
    :rtype: 2D array of float
    """
    if hom_axis is None:
        x = x / np.sqrt(np.mean(x ** 2))
    elif hom_axis == 'x':
        x = x / np.sqrt(np.mean(x ** 2, axis=1))
    elif hom_axis == 'y':
        x = x / np.sqrt(np.mean(x ** 2, axis=0))
    else:
        sys.exit('Invalid homogenity axis.')
    return x


def window(vfield, x_center_index, y_center_index, dist):
    """
    Defines a window around (x; y) coordinates

    :param vfield: fullsize velocity field
    :type vfield: 2D array of float
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
    x_index, y_index = np.meshgrid(vfield.x_coordinate_matrix[int(x1):int(x2)],
                                   vfield.y_coordinate_matrix[int(y1):int(y2)],
                                   indexing='xy')
    u_data = vfield.u_velocity_matrix[int(y1):int(y2), int(x1):int(x2)]
    v_data = vfield.v_velocity_matrix[int(y1):int(y2), int(x1):int(x2)]
    return x_index, y_index, u_data, v_data


def find_peaks(data, threshold, box_size):
    """
    Find local peaks in an image that are above above a specified
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
    :type data: 2D array of float
    :type threshold: float
    :type box_size: int

    :returns: An array containing the x and y pixel location of the peaks and their values.
    :rtype: list
    """

    if np.all(data == data.flat[0]):
        return []

    data_max = sp.ndimage.maximum_filter(data, size=box_size,
                                      mode='constant', cval=0.0)

    peak_goodmask = (data == data_max)  # good pixels are True

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[y_peaks, x_peaks]
    peaks = (y_peaks, x_peaks, peak_values)
    return peaks


def direction_rotation(vorticity, peaks):
    """
    Identify the direction of the vortices rotation using the vorticity.

    :param vorticity: 2D array with the computed vorticity
    :param peaks: list of the detected peaks
    :type vorticity: 2D array of float
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
    vortices_counterclockwise = (vortices_counterclockwise_x, vortices_counterclockwise_y,
                                 vortices_counterclockwise_cpt)
    vortices_clockwise = np.asarray(vortices_clockwise)
    vortices_counterclockwise = np.asarray(vortices_counterclockwise)
    return vortices_clockwise, vortices_counterclockwise


def correlation_coef(u_data, v_data, u, v):
    """Calculates the correlation coefficient between two 2D arrays

    :param u_data: velocity u from the data at the proposed window
    :param v_data: velocity v from the data at the proposed window
    :param u: velocity u from the calculated model
    :param v: velocity v from the calculated model
    :type u_data: float
    :type v_data: float
    :type u: float
    :type v: float
    :returns: correlation
    :rtype: float
    """
    u_data = u_data.ravel()
    v_data = v_data.ravel()
    u = u.ravel()
    v = v.ravel()

    prod_piv_mod = np.mean(u_data * u + v_data * v)
    prod_piv = np.mean(u * u + v * v)
    prod_mod = np.mean(u_data * u_data + v_data * v_data)
    correlation = prod_piv_mod / (max(prod_piv, prod_mod))

    return correlation


def velocity_model(core_radius, gamma, x_real, y_real, u_advection, v_advection, x, y):
    """Generates the Lamb-Oseen vortex velocity array

    :param core_radius: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x_real: relative x position of the vortex center
    :param y_real: relative y position of the vortex center
    :param u_advection: u advective velocity at the center
    :param v_advection: v advective velocity at the center
    :param x:
    :param y:
    :type core_radius: float
    :type gamma: float
    :type x_real: float
    :type y_real: float
    :type u_advection: float
    :type v_advection: float
    :type x: float
    :type y: float
    :returns: velx, vely
    :rtype: float
    """
    r = np.hypot(x - x_real, y - y_real)
    vel = (gamma / (2 * np.pi * r)) * (1 - np.exp(-(r ** 2) / core_radius ** 2))
    vel = np.nan_to_num(vel)
    velx = u_advection - vel * (y - y_real) / r
    vely = v_advection + vel * (x - x_real) / r
    velx = np.nan_to_num(velx)
    vely = np.nan_to_num(vely)
    # print(core_radius, gamma, x_real, y_real, u_advection, v_advection, x, y)
    return velx, vely


def get_vortices(vfield, peaks, vorticity, rmax, correlation_treshold):
    """
    General routine to check if the detected vortex is a real vortex

    :param vfield: data from the input file
    :param peaks: list of vortices
    :param vorticity: calculated field
    :param rmax: maximum radius (adapt it to your data domain)
    :param correlation_treshold: threshold to detect a vortex (default is 0.75)
    :type vfield: class VelocityField
    :type peaks: list
    :type vorticity: array
    :type rmax: float
    :type correlation_treshold: float
    :returns: list of detected vortices
    :rtype: list
    """

    vortices = list()
    cpt_accepted = 0
    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step
    for i in range(len(peaks[0])):
        x_center_index = peaks[1][i]
        y_center_index = peaks[0][i]
        print(i, 'Processing detected swirling at (x, y)', x_center_index, y_center_index)
        if rmax == 0.0:
            core_radius = 2 * np.hypot(dx, dy)
        else:
            core_radius = rmax  # guess on the starting vortex radius
        gamma = vorticity[y_center_index, x_center_index] * np.pi * core_radius ** 2

        vortices_parameters = full_fit(core_radius, gamma, vfield, x_center_index, y_center_index)
        if vortices_parameters[6] < 2:
            correlation_value = 0
        else:
            x_index, y_index, u_data, v_data = window(vfield, round(vortices_parameters[2] / dx, 0),
                                                            round(vortices_parameters[3] / dy, 0),
                                                            vortices_parameters[6])
            u_model, v_model = velocity_model(vortices_parameters[0], vortices_parameters[1], vortices_parameters[2],
                                              vortices_parameters[3],
                                              vortices_parameters[4], vortices_parameters[5], x_index, y_index)
            correlation_value = correlation_coef(u_data - vortices_parameters[4], v_data - vortices_parameters[5],
                                                 u_model - vortices_parameters[4], v_model - vortices_parameters[5])
        if correlation_value > correlation_treshold:
            print('Accepted! Correlation = {:1.2f} (vortex #{:2d})'.format(correlation_value, cpt_accepted))
            u_theta = (vortices_parameters[1] / (2 * np.pi * vortices_parameters[0])) * (
                    1 - np.exp(-1))  # compute the tangential velocity at critical radius
            vortices.append(
                [vortices_parameters[0], vortices_parameters[1], vortices_parameters[2], vortices_parameters[3],
                 vortices_parameters[4],
                 vortices_parameters[5], vortices_parameters[6], correlation_value, u_theta])
            cpt_accepted += 1
    return vortices


def full_fit(core_radius, gamma, vfield, x_center_index, y_center_index):
    """Full fitting procedure

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
    :returns: fitted[i], dist
    :rtype: list
    """

    fitted = [[], [], [], [], [], []]
    fitted[0] = core_radius
    fitted[1] = gamma
    fitted[2] = vfield.x_coordinate_matrix[x_center_index]
    fitted[3] = vfield.y_coordinate_matrix[y_center_index]
    dx = vfield.x_coordinate_step
    dy = vfield.y_coordinate_step
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
        dist = int(round(fitted[0] / np.hypot(dx, dy), 0)) + 1
        if fitted[0] < 2 * np.hypot(dx,dy):
            break
        fitted[4] = vfield.u_velocity_matrix[y_center_index, x_center_index]  # u_advection
        fitted[5] = vfield.v_velocity_matrix[y_center_index, x_center_index]  # v_advection
        x_index, y_index, u_data, v_data = window(vfield, x_center_index, y_center_index, dist)


        fitted = fit(fitted[0], fitted[1], x_index, y_index, fitted[2], fitted[3],
                     u_data, v_data, fitted[4], fitted[5], i)
        if i > 0:
            # break if radius variation is less than 10% and accepts
            if abs(fitted[0] / r1 - 1) < 0.1:
                if (abs((fitted[2] / x1 - 1)) < 0.1) or (abs((fitted[3] / y1 - 1)) < 0.1):
                    break
            # break if x or y position is out of the window and discards
            if (abs((fitted[2] - x1)) > dist * dx) or (abs((fitted[3] - y1)) > dist * dy):
                dist = 0
                break
    return fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], dist


def fit(core_radius, gamma, x, y, x_real, y_real, u_data, v_data, u_advection, v_advection, i):
    """
    Fitting  of the Lamb-Oseen Vortex

    :param core_radius: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x: x position
    :param y: y position
    :param x_real: x position of the vortex center
    :param y_real: y position of the vortex center
    :type core_radius: float
    :type gamma: float
    :type x: float
    :type y: float
    :type x_real: float
    :type y_real: float
    :returns: sol.x
    :rtype: float
    """

    x = x.ravel()
    y = y.ravel()
    u_data = u_data.ravel()
    v_data = v_data.ravel()
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    def lamb_oseen_model(fitted):
        """
        Lamb-Oseen velocity model used for the nonlinear fitting
        """
        r = np.hypot(x - fitted[2], y - fitted[3])
        expr2 = np.exp(-r ** 2 / fitted[0] ** 2)
        z = fitted[1] / (2 * np.pi * r) * (1 - expr2)
        z = np.nan_to_num(z)
        zx = fitted[4] - z * (y - fitted[3]) / r - u_data
        zy = fitted[5] + z * (x - fitted[2]) / r - v_data
        zx = np.nan_to_num(zx)
        zy = np.nan_to_num(zy)
        zt = np.append(zx, zy)
        return zt

    if i > 0:
        m = 1.0
    else:
        m = 4.0

    epsilon = 0.001
    bnds = ([0, gamma - abs(gamma) * m / 2 - epsilon, x_real - m * dx - epsilon, y_real - m * dy - epsilon,
             u_advection - abs(u_advection) - epsilon, v_advection - abs(v_advection) - epsilon],
            [core_radius + core_radius * m, gamma + abs(gamma) * m / 2 + epsilon, x_real + m * dx + epsilon,
             y_real + m * dy + epsilon, u_advection + abs(u_advection) + epsilon, v_advection + abs(v_advection) + epsilon])

    sol = opt.least_squares(lamb_oseen_model, [core_radius, gamma, x_real, y_real, u_advection, v_advection],
                                    method='trf', bounds=bnds)

    return sol.x


def plot_fields(vfield, detection_field):
    """
    Plot fields: display the (u,v,w) fields and the vorticity field.

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :returns: popup
    :rtype: image
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)  # , sharex='col', sharey='row')
    ax1.imshow(vfield.u_velocity_matrix, cmap='seismic', origin="lower")
    ax1.set_title('Velocity u (velocity_s)')

    ax2.imshow(vfield.v_velocity_matrix, cmap='seismic', origin="lower")
    ax2.set_title('Velocity v (velocity_n)')

    try:
        ax3.imshow(vfield.w_velocity_matrix, cmap='seismic', origin="lower")
        ax3.set_title('Velocity w (velocity_z)')
    except:
        print('No w velocity')
    ax4.set_title('Vorticity')
    ax4.imshow(detection_field, origin="lower", cmap='seismic')
    plt.tight_layout()

    plt.show()


def plot_detect(vortices_counterclockwise, vortices_clockwise, detection_field, *args):
    """
    Plot detect: display the location and rotation of the vortices

    :param vortices_counterclockwise: vortices spinning counterclockwise
    :type vortices_counterclockwise: list of vortices
    :param vortices_clockwise: vortices spinning clockwise
    :type vortices_clockwise: list of vortices
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :returns: popup
    :rtype: image
    """
    plt.subplot()
    if (args[0] == True):
        detection_field = detection_field.T  # transpose the detection field
        plt.scatter(vortices_counterclockwise[0], vortices_counterclockwise[1], edgecolor='G', facecolor='G',
                    label='left')
        plt.scatter(vortices_clockwise[0], vortices_clockwise[1], edgecolor='Y', facecolor='Y', label='right')
    else:
        plt.scatter(vortices_counterclockwise[1], vortices_counterclockwise[0], edgecolor='G', facecolor='G',
                    label='left')
        plt.scatter(vortices_clockwise[1], vortices_clockwise[0], edgecolor='Y', facecolor='Y', label='right')

    plt.title('Detected possible vortices')
    # plt.contourf(field, cmap="Greys_r")

    plt.imshow(detection_field, origin='lower', cmap="Greys_r")
    plt.xlabel('x')
    plt.ylabel('y')
    # plt.legend()
    # plt.imshow(field, cmap="Greys_r",origin="lower")
    plt.tight_layout()

    plt.show()


def plot_quiver(x_index, y_index, u_data, v_data, detection_field):
    """
    Plot quiver: display a specific (x,y) location with vector fields.

    :param x_index: contains spatial mesh (x direction)
    :type x_index: array of float
    :param y_index: contains spatial mesh (y direction)
    :type y_index: array of float
    :param u_data: contains velocity data (u component)
    :type u_data: 2D array of float
    :param v_data: contains velocity data (v component)
    :type v_data: 2D array of float
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :returns: popup
    :rtype: image
    """

    plt.figure()
    # plt.title('Velocity vectors centered at max swirling strength')
    plt.contourf(detection_field,
                 extent=[x_index[0][0], x_index[0][-1], y_index[0][0], y_index[-1][0]])
    s = 1  # sampling factor, can be modified
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s])
    plt.show()


def plot_fit(x_index, y_index, u_data, v_data, u_model, v_model,
             xc, yc, core_radius, gamma, u_advection, v_advection, correlation_value,
             cpt_vortex, subtract_advection_field, output_dir, time_step):
    """
    Plot fit

    :param x_index: contains spatial mesh (x direction)
    :type x_index: array of float
    :param y_index: contains spatial mesh (y direction)
    :type y_index: array of float
    :param u_data: contains velocity data (u component)
    :type u_data: 2D array of float
    :param v_data: contains velocity data (v component)
    :type v_data: 2D array of float
    :param u_model: contains velocity data (u component)
    :type u_model: 2D array of float
    :param v_model: contains velocity data (v component)
    :type v_model: 2D array of float
    :param xc: x coordinate of the vortex center
    :type xc: float
    :param yc: y coordinate of the vortex center
    :type yc: float
    :param core_radius: dimension of the vortex core radius
    :type core_radius: float
    :param gamma: circulation of the vortex
    :type gamma: float
    :param u_advection: contains velocity data (u component)
    :type u_advection: 2D array of float
    :param v_advection: contains velocity data (v component)
    :type v_advection: 2D array of float
    :param cpt_vortex: current n° of the vortex
    :type cpt_vortex: int
    :param subtract_advection_field: if True, the advection field (u_advection, v_advection) is subtracted
    :type subtract_advection_field:  bool
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param correlation_value: correlation between the vortex and a Lamb-Oseen model
    :type correlation_value: float
    :param time_step: current time_step
    :type time_step: int
    :returns: image file
    :rtype: png
    """
    plt.figure()
    s = 1  # sampling factor
    if (x_index.size > 400):
        s = 2
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_data[::s, ::s], v_data[::s, ::s],
               color='r', label='data')
    plt.quiver(x_index[::s, ::s], y_index[::s, ::s], u_model[::s, ::s], v_model[::s, ::s],
               color='b', label='model', alpha=0.5)
    circle1 = plt.Circle((xc, yc), core_radius, color='k', alpha=0.05)
    plt.gca().add_artist(circle1)
    plt.gca().scatter([xc], [yc], marker='+', color='k', s=100)
    plt.legend()
    plt.grid()
    plt.gca().set_aspect('equal', adjustable='box')
    #    plt.axes().set_aspect('equal') #deprecated
    plt.xlabel('x')
    plt.ylabel('y')

    plt.title(r'r=%s $\Gamma$=%s u=%s v=%s C=%s' % (
        round(core_radius, 2), round(gamma, 2), round(u_advection, 2), round(v_advection, 2),
        round(correlation_value, 2)))
    if subtract_advection_field == False:
        plt.savefig(output_dir + '/vortex%i_%s.png' % (cpt_vortex, 'initial_vfield'), format='png')
    else:
        plt.savefig(output_dir + '/vortex%i_%s.png' % (cpt_vortex, 'advection_field_subtracted'), format='png')
    plt.close('all')


def plot_accepted(vfield, vortices_list, detection_field, output_dir, time_step):
    """
    Plot accepted: display the accepted vortices, with respect to the different criterion
    (correlation threshold, box size ...)

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param detection_field: detection field (vorticity ...)
    :type detection_field: 2D array of float
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
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
    for i, line in enumerate(vortices_list):
        if vortices_list[i][1] > 0:
            orient = 'Y'
        else:
            orient = 'Y'
        plt.figure(1)
        circle1 = plt.Circle((line[2], line[3]), line[0],
                             edgecolor='yellow', facecolor='none', gid='vortex%i' % i)
        plt.gca().add_artist(circle1)
        plt.figure(2)
        circle1 = plt.Circle((line[2] / dx, line[3] / dy), line[0] / np.hypot(dx, dy),
                             edgecolor='yellow', facecolor='none', gid='vortex%i' % i)
        plt.gca().add_artist(circle1)

    ##Comparing data
    # fileIn = open('../data/dazin.dat', 'r')
    # for line in fileIn:
    # xComp = int(float(line.split()[1]))
    # yComp = int(float(line.split()[2]))
    # gammaComp = float(line.split()[3])
    # rComp = float(line.split()[4])
    # if gammaComp > 0:
    # orient = 'R'
    # else:
    # orient = 'R'
    # circle2=plt.Circle((xComp,yComp),rComp,edgecolor=orient,facecolor='none')
    # plt.gca().add_artist(circle2)

    # plt.legend()
    plt.figure(1)
    plt.tight_layout()
    plt.savefig(output_dir + '/accepted_{:01d}.svg'.format(time_step), format='svg')
    create_links(output_dir + '/accepted_{:01d}.svg'.format(time_step), vortices_list, output_dir, time_step)
    plt.figure(2)
    plt.savefig(output_dir + '/meshed_{:01d}.svg'.format(time_step),
                format='svg')  # use it to verify your coordinate system if needed !
    # plt.savefig(output_dir+'/tk_{:01d}.png'.format(time_step), format='png', transparent=True)

    # plt.show()


def plot_vortex(vfield, vortices_list, output_dir, time_step):
    """
    Plot vortex: plot a vortex and its corresponding vortex model

    :param vfield: contains spatial mesh and velocity components
    :type vfield: class VelocityField()
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :returns: file
    :rtype: image
    """
    for cpt_vortex, line in enumerate(vortices_list):
        print('r: {:.3f}'.format(line[0]),
              'gamma: {:.2f}'.format(line[1]),
              'x: {:.2f}'.format(line[2]),
              'y: {:.2f}'.format(line[3]),
              'dist: {:.2f}'.format(line[6]),
              'corr: {:.2f}'.format(line[7]),
              'vt: {:.2f}'.format(line[8]))
        dx = vfield.x_coordinate_step
        dy = vfield.y_coordinate_step
        x_index, y_index, u_data, v_data = window(vfield, round(line[2] / dx, 0), round(line[3] / dy, 0), line[6])
        u_model, v_model = velocity_model(line[0], line[1], line[2], line[3], line[4], line[5], x_index,
                                          y_index)
        correlation_value = correlation_coef(u_data, v_data, u_model, v_model)
        plot_fit(x_index, y_index, u_data, v_data, u_model, v_model, line[2], line[3],
                 line[0], line[1], line[4], line[5], correlation_value, cpt_vortex, False, output_dir, time_step)
        correlation_value = correlation_coef(u_data - line[4], v_data - line[5], u_model - line[4], v_model - line[5])
        plot_fit(x_index, y_index, u_data - line[4], v_data - line[5], u_model - line[4], v_model - line[5], line[2],
                 line[3],
                 line[0], line[1], line[4], line[5], correlation_value, cpt_vortex, True, output_dir, time_step)


def create_links(path, vortices_list, output_dir, time_step):
    """
    create links: add some links bewteen the accepted.svg file and the detected vortices

    :param path: path of the accepted.svg file
    :type path: str
    :param vortices_list: contains all the detected vortices
    :type vortices_list: list
    :param output_dir: directory where the results are written
    :type output_dir: str
    :param time_step: current time_step
    :type time_step: int
    :returns: file
    :rtype: image
    """

    file_in = open(path, "r")
    file_out = open(output_dir + "/linked_{:01d}.svg".format(time_step), "w")
    i = 0
    vortex_found = False
    for line in file_in:
        if "</g>" in line:
            if vortex_found == True:
                file_out.write(line)
                file_out.write('   </a>\n')
                vortex_found = False
            else:
                file_out.write(line)
        elif "vortex" in line:
            file_out.write('   <a href="vortex%i_advection_field_subtracted.png">\n' % i)
            file_out.write(line)
            file_out.write('   <title>Vortex %i: r = %s gamma = %s</title>\n' % (
                i, round(vortices_list[i][3], 1), round(vortices_list[i][2], 1)))
            i = i + 1
            vortex_found = True
        else:
            file_out.write(line)
