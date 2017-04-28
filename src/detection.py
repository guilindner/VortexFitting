import numpy as np
from scipy import ndimage

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
        detection threshold.  A 2D ``threshold`` must have the same
        shape as ``data``.  See `detect_threshold` for one way to create
        a ``threshold`` image.

    box_size : scalar or tuple, optional
        The size of the local region to search for peaks at every point
        in ``data``.  If ``box_size`` is a scalar, then the region shape
        will be ``(box_size, box_size)``.  Either ``box_size`` or
        ``footprint`` must be defined.  If they are both defined, then
        ``footprint`` overrides ``box_size``.

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
    y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[y_peaks, x_peaks]

    peaks = (x_peaks, y_peaks, peak_values)

    return peaks

def direction_rotation(vorticity,peaks):
    """ Identify the direction of the vortices rotation
    using the vorticity.
    """
    clockwise = []
    clockwise_x, clockwise_y, clockwise_i = [],[],[]
    counterclockwise = []
    counterclockwise_x, counterclockwise_y, counterclockwise_i = [],[],[]
    for i in range(len(peaks[0])):
        if vorticity[peaks[1][i],peaks[0][i]] > 0.0:
            clockwise_x.append(peaks[0][i])
            clockwise_y.append(peaks[1][i])
            clockwise_i.append(peaks[2][i])
        else:
            counterclockwise_x.append(peaks[0][i])
            counterclockwise_y.append(peaks[1][i])
            counterclockwise_i.append(peaks[2][i])
    clockwise = (clockwise_x, clockwise_y, clockwise_i)
    counterclockwise = (counterclockwise_x, counterclockwise_y, counterclockwise_i)
    clockwise = np.asarray(clockwise)
    counterclockwise = np.asarray(counterclockwise)
    return clockwise, counterclockwise
