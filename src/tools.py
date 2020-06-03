import sys
import numpy as np
from scipy import ndimage

np.seterr(divide='ignore', invalid='ignore')


def get_fluc(x, mean, hom_axis):
    """
    Used when you have a advective velocity along one axis

    :param x: velocity field
    :type x: 2D array of float
    :param mean: advective velocity to subtract
    :type mean: float
    :param hom_axis: False, 'x', or 'y'. The axis which the mean is subtracted
    :type hom_axis: str
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
    """
    if hom_axis is None:
        x = x/np.sqrt(np.mean(x**2))
    elif hom_axis == 'x':
        x = x/np.sqrt(np.mean(x**2, axis=1))
    elif hom_axis == 'y':
        x = x/np.sqrt(np.mean(x**2, axis=0))
    else:
        sys.exit("Invalid homogenity axis.")
    return x

def window(vfield,x_center_index,y_center_index,dist):
    """
    Defines a window around (x;y) coordinates

    :param vfield: fullsize velocity field
    :type vfield: 2D array of float
    :param x_center_index: box center index (x)
    :type x_center_index: int
    :param y_center_index: box center index (y)
    :type y_center_index: int
    :param dist: ??    
    :param dist: int
    :returns: cropped arrays for x, y, u and v 
    :rtype: 2D arrays of floats

    """
    if (x_center_index-dist > 0):
        x1 = x_center_index -dist
    else:
        x1 = 0
    if (y_center_index-dist > 0):
        y1 = y_center_index -dist
    else:
        y1 = 0
    if (x_center_index+dist <= vfield.u.shape[1]):
        x2 = x_center_index+dist
    else:
        x2 = vfield.u.shape[0]
    if (y_center_index+dist <= vfield.v.shape[0]):
        y2 = y_center_index+dist
    else:
        y2 = vfield.v.shape[0]
    x_index, y_index = np.meshgrid(vfield.dx[int(x1):int(x2)],
                       vfield.dy[int(y1):int(y2)],indexing='xy')
    u_data = vfield.u[int(y1):int(y2),int(x1):int(x2)]
    v_data = vfield.v[int(y1):int(y2),int(x1):int(x2)]
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

    data_max = ndimage.maximum_filter(data, size=box_size,
                                          mode='constant', cval=0.0)

    peak_goodmask = (data == data_max)    # good pixels are True

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[y_peaks, x_peaks]
    peaks = (y_peaks, x_peaks, peak_values)
    return peaks

def direction_rotation(vorticity,peaks):
    """ Identify the direction of the vortices rotation using the vorticity.
    :param vorticity: 2D array with the computed vorticity
    :param peaks: list of the detected peaks

    :returns: An array containing the direction of rotation for each vortex
    :rtype: list
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
