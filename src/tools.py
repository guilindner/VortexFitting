import numpy as np
from classes import VelocityField
from scipy import ndimage

np.seterr(divide='ignore', invalid='ignore')

def sub_mean(x, hom_axis):
    """
    Used when you have a convective velocity along one axis
    """
    mean = np.mean(x, axis=hom_axis)
    if hom_axis == 0:
        x = x - mean
    else:
        x = x - mean[:, None]
    return x

def normalize(x, hom_axis):
    mean = np.mean(x**2, axis=hom_axis)
    x = x/np.sqrt(mean)
    return x

def window(a,x_center_index,y_center_index,dist):
    if (x_center_index-dist > 0):
        x1 = x_center_index -dist
    else:
        x1 = 0
    if (y_center_index-dist > 0):
        y1 = y_center_index -dist
    else:
        y1 = 0
    if (x_center_index+dist <= a.u.shape[1]):
        x2 = x_center_index+dist
    else:
        x2 = a.u.shape[1]
    if (y_center_index+dist <= a.v.shape[0]):
        y2 = y_center_index+dist
    else:
        y2 = a.v.shape[0]
    x_index, y_index = np.meshgrid(a.dx[int(x1):int(x2)],
                       a.dy[int(y1):int(y2)],indexing='xy')
    u_data = a.u[int(y1):int(y2),int(x1):int(x2)]
    v_data = a.v[int(y1):int(y2),int(x1):int(x2)]
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

    :returns: An array containing the x and y pixel location of the peaks and their values.
    :rtype: list
    """

    if np.all(data == data.flat[0]):
        return []

    data_max = ndimage.maximum_filter(data, size=box_size,
                                          mode='constant', cval=0.0)

    peak_goodmask = (data == data_max)    # good pixels are True

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[y_peaks, x_peaks]
    peaks = (y_peaks, x_peaks, peak_values)
    return peaks

def direction_rotation(vorticity,peaks):
    """ Identify the direction of the vortices rotation
    using the vorticity.
    """
    dirR = []
    dirR_x, dirR_y, dirR_i = [],[],[]
    dirL = []
    dirL_x, dirL_y, dirL_i = [],[],[]
    for i in range(len(peaks[0])):
        if vorticity[peaks[0][i],peaks[1][i]] > 0.0:
            dirR_x.append(peaks[0][i])
            dirR_y.append(peaks[1][i])
            dirR_i.append(peaks[2][i])
        else:
            dirL_x.append(peaks[0][i])
            dirL_y.append(peaks[1][i])
            dirL_i.append(peaks[2][i])
    dirR = (dirR_x, dirR_y, dirR_i)
    dirL = (dirL_x, dirL_y, dirL_i)
    dirR = np.asarray(dirR)
    dirL = np.asarray(dirL)
    return dirR, dirL
