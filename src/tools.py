import numpy as np
from classes import VelocityField
from scipy import ndimage

np.seterr(divide='ignore', invalid='ignore')

def sub_mean(x, hom_axis):
    mean = np.mean(x, axis=hom_axis)
    if (hom_axis == 0):
        x = x - mean
    else:
        x = x - mean[:,None]
    return x
      
def normalize(x, hom_axis):
    mean = np.mean(x**2, axis=hom_axis)
    x = x/np.sqrt(mean)
    return x

def window(a,xCenter,yCenter,dist):
    if (xCenter-dist >= 0):
        x1 = xCenter -dist
    else:
        x1 = 0
    if (yCenter-dist >= 0):
        y1 = yCenter -dist
    else:
        y1 = 0
    if (xCenter+dist <= len(a.u[0]-1)):
        x2 = xCenter+dist
    else:
        x2 = len(a.u[0]-1)
    if (yCenter+dist <= len(a.v[0]-1)):
        y2 = yCenter+dist
    else:
        y2 = len(a.v[0]-1)
    X, Y = np.meshgrid(a.dx[x1:x2],
                       a.dy[y1:y2],indexing='xy')
    Uw = a.u[x1:x2,y1:y2]
    Vw = a.v[x1:x2,y1:y2]
    return X, Y, Uw, Vw

def find_peaks(data, threshold, box_size):
    """
    Find local peaks in an image that are above above a specified
    threshold value.

    Peaks are the maxima above the "threshold" within a local region.
    The regions are defined by the "box_size" parameters.
    "box_size" defines the local region around each pixel
    as a square box.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    threshold : float or array-like
        The data value or pixel-wise data values to be used for the
        detection threshold.  A 2D "threshold" must have the same
        shape as "data".

    box_size : scalar or tuple, optional
        The size of the local region to search for peaks at every point

    Returns
    -------
    output :
        An array containing the x and y pixel location of the peaks and
        their values.
    """

    if np.all(data == data.flat[0]):
        return []

    data_max = ndimage.maximum_filter(data, size=box_size,
                                          mode='constant', cval=0.0)

    peak_goodmask = (data == data_max)    # good pixels are True

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    x_peaks, y_peaks = peak_goodmask.nonzero()
    peak_values = data[x_peaks, y_peaks]

    peaks = (x_peaks, y_peaks, peak_values)

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
