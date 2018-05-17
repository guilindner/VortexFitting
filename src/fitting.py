import numpy as np
from scipy.signal import correlate2d
from scipy import optimize
from scipy.stats import pearsonr

import tools
import plot

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
    :returns: corr
    :rtype: float
    """
    u_data = u_data.ravel()
    v_data = v_data.ravel()
    u = u.ravel()
    v = v.ravel()
    N = u_data.size

    prod_PIV_mod = np.mean(u_data*u + v_data*v)
    prod_PIV     = np.mean(u*u      + v*v)
    prod_mod     = np.mean(u_data*u_data + v_data*v_data)
    corr = prod_PIV_mod/(max(prod_PIV,prod_mod))

    return corr

def velocity_model(coreR, gamma, x_real, y_real, u_conv, v_conv, x, y):
    """Generates the Lamb-Oseen vortex velocity array

    :param coreR: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x_real: relative x position of the vortex center
    :param y_real: relative y position of the vortex center
    :param u_conv: u convective velocity at the center
    :param v_conv: v convective velocity at the center
    :type coreR: float
    :type gamma: float
    :type x_real: float
    :type y_real: float
    :type u_conv: float
    :type v_conv: float
    :returns: velx, vely
    :rtype: float
    """
    r = np.hypot(x-x_real, y-y_real)
    vel = (gamma/(2 * np.pi * r)) * (1 - np.exp(-(r**2)/(coreR)**2))
    vel = np.nan_to_num(vel)
    velx = u_conv -vel*(y-y_real)/r
    vely = v_conv +vel*(x-x_real)/r
    velx = np.nan_to_num(velx)
    vely = np.nan_to_num(vely)
    return velx, vely

def get_vortices(a, peaks, vorticity,rmax):
    """General routine to check if the detected vortex is a real vortex

    :param a: data from the input file
    :param peaks: list of vortices
    :param vorticity: calculated field
    :type a: int, float, ...
    :type peaks: list
    :type vorticity: array
    :returns: vortices
    :rtype: list
    """
    b = list()
    vortices = list()
    j = 0
    dx = a.dx[5]-a.dx[4] #ugly
    dy = a.dy[5]-a.dy[4]
    for i in range(len(peaks[0])):
        x_center_index = peaks[1][i]
        y_center_index = peaks[0][i]
        print(i, " Processing detected swirling at (x, y)", x_center_index, y_center_index)
        #coreR = 2*(a.dx[5]-a.dx[4]) #ugly change someday
        coreR = rmax #guess on the starting vortex radius
        gamma = vorticity[y_center_index, x_center_index]*np.pi*coreR**2
        #print("vorticity",vorticity[y_center_index, x_center_index])
        b = full_fit(coreR, gamma, a, x_center_index, y_center_index)
        if b[6] < 2:
            corr = 0
        else:
            x_index, y_index, u_data, v_data = tools.window(a, round(b[2]/dx, 0), round(b[3]/dy, 0), b[6])
            u_model, v_model = velocity_model(b[0], b[1], b[2], b[3], b[4], b[5], x_index, y_index)
            corr = correlation_coef(u_data-b[4], v_data-b[5], u_model-b[4], v_model-b[5])
        if corr > 0.55: #if the vortex is too big, its better to decrease this value
            print("Accepted! corr = %s (vortex %s)" %(corr, j))
            velT = (b[1]/(2 * np.pi * b[0])) * (1 - np.exp(-1))
            vortices.append([b[0], b[1], b[2], b[3], b[4], b[5], b[6], corr, velT])
            j += 1
    return vortices

def full_fit(coreR, gamma, a, x_center_index, y_center_index):
    """Full fitting procedure

    :param coreR: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param a: data from the input file
    :param x_center_index: x index of the vortex center
    :param y_center_index: y index of the vortex center
    :type coreR: float
    :type gamma: float
    :type a: class
    :type x_center_index: int
    :type y_center_index: int
    :returns: fitted[i], dist
    :rtype: list
    """

    fitted = [[], [], [], [], [], []]
    fitted[0] = coreR
    fitted[1] = gamma
    fitted[2] = a.dx[x_center_index]
    fitted[3] = a.dy[y_center_index]
    dx = a.dx[5]-a.dx[4] #ugly
    dy = a.dy[5]-a.dy[4]
    corr = 0.0
    for i in range(5):
        x_center_index = int(round(fitted[2]/dx))
        y_center_index = int(round(fitted[3]/dy))
        if x_center_index >= a.u.shape[1]:
            x_center_index = a.u.shape[1]-1
        if x_center_index <= 2:
            x_center_index = 3
        if y_center_index >= a.v.shape[0]:
            y_center_index = a.v.shape[0]-1
        r1 = fitted[0]
        x1 = fitted[2]
        y1 = fitted[3]
        dist = int(round(fitted[0]/dx, 0)) + 1
        if fitted[0] < 2*dx:
            break
        fitted[4] = a.u[y_center_index, x_center_index] #u_conv
        fitted[5] = a.v[y_center_index, x_center_index] #v_conv
        x_index, y_index, u_data, v_data = tools.window(a, x_center_index, y_center_index, dist)
        fitted = fit(fitted[0], fitted[1], x_index, y_index, fitted[2], fitted[3],
                     u_data, v_data, fitted[4], fitted[5], i)
        if i > 0:
            # break if radius variation is less than 10% and accepts
            if abs(fitted[0]/r1 -1) < 0.1:
                if (abs((fitted[2]/x1 -1)) < 0.1) or (abs((fitted[3]/y1 -1)) < 0.1):
                    break
            # break if x or y position is out of the window and discards
            if (abs((fitted[2]-x1)) > dist*dx) or (abs((fitted[3]-y1)) > dist*dy):
                dist = 0
                break
    return fitted[0], fitted[1], fitted[2], fitted[3], fitted[4], fitted[5], dist
def fit(coreR, gamma, x, y, x_real, y_real, u_data, v_data, u_conv, v_conv, i):
    """
    Fitting  of the Lamb-Oseen Vortex

    :param coreR: core radius of the vortex
    :param gamma: circulation contained in the vortex
    :param x: x position
    :param y: y position
    :param x_real: x position of the vortex center
    :param y_real: y position of the vortex center
    :type coreR: float
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
    dx = x[1]-x[0]
    dy = dx
    def fun(fitted):
        """
        Lamb-Ossen velocity model used for the nonlinear fitting
        """
        r = np.hypot(x-fitted[2], y-fitted[3])
        expr2 = np.exp(-r**2/fitted[0]**2)
        z = fitted[1]/(2*np.pi*r) * (1 - expr2)
        z = np.nan_to_num(z)
        zx = fitted[4] -z*(y-fitted[3])/r -u_data
        zy = fitted[5] +z*(x-fitted[2])/r -v_data
        zx = np.nan_to_num(zx)
        zy = np.nan_to_num(zy)
        zt = np.append(zx, zy)
        return zt
    if i > 0:
        m = 1.0
    else:
        m = 4.0
    bnds = ([0.0, gamma-abs(gamma)*m/2, x_real-m*dx, y_real-m*dy,
             u_conv-abs(u_conv), v_conv-abs(v_conv)],
            [coreR+coreR*m, gamma+abs(gamma)*m/2, x_real+m*dx, y_real+m*dy,
             u_conv+abs(u_conv), v_conv+abs(v_conv)])
    sol = optimize.least_squares(fun, [coreR, gamma, x_real, y_real, u_conv,
                                       v_conv], bounds=bnds)
    return sol.x
