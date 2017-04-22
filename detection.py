import numpy as np

def find_peaks(data, threshold, box_size=3, footprint=None, mask=None,
               border_width=None, npeaks=np.inf, error=None,
               wcs=None):
    """
    Find local peaks in an image that are above above a specified
    threshold value.

    Peaks are the maxima above the ``threshold`` within a local region.
    The regions are defined by either the ``box_size`` or ``footprint``
    parameters.  ``box_size`` defines the local region around each pixel
    as a square box.  ``footprint`` is a boolean array where `True`
    values specify the region shape.

    If multiple pixels within a local region have identical intensities,
    then the coordinates of all such pixels are returned.  Otherwise,
    there will be only one peak pixel per local region.  Thus, the
    defined region effectively imposes a minimum separation between
    peaks (unless there are identical peaks within the region).

    When using subpixel precision (``subpixel=True``), then a cutout of
    the specified ``box_size`` or ``footprint`` will be taken centered
    on each peak and fit with a 2D Gaussian (plus a constant).  In this
    case, the fitted local centroid and peak value (the Gaussian
    amplitude plus the background constant) will also be returned in the
    output table.

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

    footprint : `~numpy.ndarray` of bools, optional
        A boolean array where `True` values describe the local footprint
        region within which to search for peaks at every point in
        ``data``.  ``box_size=(n, m)`` is equivalent to
        ``footprint=np.ones((n, m))``.  Either ``box_size`` or
        ``footprint`` must be defined.  If they are both defined, then
        ``footprint`` overrides ``box_size``.

    mask : array_like, bool, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.

    border_width : bool, optional
        The width in pixels to exclude around the border of the
        ``data``.

    npeaks : int, optional
        The maximum number of peaks to return.  When the number of
        detected peaks exceeds ``npeaks``, the peaks with the highest
        peak intensities will be returned.

    subpixel : bool, optional
        If `True`, then a cutout of the specified ``box_size`` or
        ``footprint`` will be taken centered on each peak and fit with a
        2D Gaussian (plus a constant).  In this case, the fitted local
        centroid and peak value (the Gaussian amplitude plus the
        background constant) will also be returned in the output table.

    error : array_like, optional
        The 2D array of the 1-sigma errors of the input ``data``.
        ``error`` is used only to weight the 2D Gaussian fit performed
        when ``subpixel=True``.

    wcs : `~astropy.wcs.WCS`
        The WCS transformation to use to convert from pixel coordinates
        to ICRS world coordinates.  If `None`, then the world
        coordinates will not be returned in the output
        `~astropy.table.Table`.

    Returns
    -------
    output : `~astropy.table.Table`
        A table containing the x and y pixel location of the peaks and
        their values.  If ``subpixel=True``, then the table will also
        contain the local centroid and fitted peak value.
    """

    from scipy import ndimage

    if np.all(data == data.flat[0]):
        return []

    if footprint is not None:
        data_max = ndimage.maximum_filter(data, footprint=footprint,
                                          mode='constant', cval=0.0)
    else:
        data_max = ndimage.maximum_filter(data, size=box_size,
                                          mode='constant', cval=0.0)

    peak_goodmask = (data == data_max)    # good pixels are True

    if mask is not None:
        mask = np.asanyarray(mask)
        if data.shape != mask.shape:
            raise ValueError('data and mask must have the same shape')
        peak_goodmask = np.logical_and(peak_goodmask, ~mask)

    if border_width is not None:
        for i in range(peak_goodmask.ndim):
            peak_goodmask = peak_goodmask.swapaxes(0, i)
            peak_goodmask[:border_width] = False
            peak_goodmask[-border_width:] = False
            peak_goodmask = peak_goodmask.swapaxes(0, i)

    peak_goodmask = np.logical_and(peak_goodmask, (data > threshold))
    y_peaks, x_peaks = peak_goodmask.nonzero()
    peak_values = data[y_peaks, x_peaks]

    if len(x_peaks) > npeaks:
        idx = np.argsort(peak_values)[::-1][:npeaks]
        x_peaks = x_peaks[idx]
        y_peaks = y_peaks[idx]
        peak_values = peak_values[idx]

    columns = (x_peaks, y_peaks, peak_values)
    names = ('x_peak', 'y_peak', 'peak_value')

#    table = Table(columns, names=names)

    if wcs is not None:
        icrs_ra_peak, icrs_dec_peak = pixel_to_icrs_coords(x_peaks, y_peaks,
                                                           wcs)
        table.add_column(Column(icrs_ra_peak, name='icrs_ra_peak'), index=2)
        table.add_column(Column(icrs_dec_peak, name='icrs_dec_peak'), index=3)

        if subpixel:
            icrs_ra_centroid, icrs_dec_centroid = pixel_to_icrs_coords(
                x_centroid, y_centroid, wcs)
            idx = table.colnames.index('y_centroid')
            table.add_column(Column(icrs_ra_centroid,
                                    name='icrs_ra_centroid'), index=idx+1)
            table.add_column(Column(icrs_dec_centroid,
                                    name='icrs_dec_centroid'), index=idx+2)

    return columns
