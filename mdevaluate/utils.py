"""
Collection of utility functions.
"""
import hashlib
import functools
from types import FunctionType

import numpy as np
import numba

from scipy.interpolate import interp1d


def hash_partial(partial):
    """
    Hashes `functools.partail` objects by their function and arguments.
    Keyword argmuents are sorted by name, to preserve order between python sessions.
    """
    hashes = [hash_anything(partial.func)]
    for arg in partial.args:
        hashes.append(hash_anything(arg))
    keys = list(partial.keywords.keys())
    keys.sort()
    for key in keys:
        hashes.append(hash_anything(partial.keywords[key]))
    return merge_hashes(*hashes)


def hash_code(func):
    """
    Hash an object with a __code__ attribute, e.g. a function.
    """
    hashes = [hash_anything(func.__code__.co_code)]
    if func.__closure__ is not None:
        for cell in func.__closure__:
            hashes.append(hash_anything(cell.cell_contents))
    return merge_hashes(*hashes)


def hash_anything(arg):
    """Return a md5 hash value for the current state of any argument."""
    if isinstance(arg, bytes):
        bstr = arg
    elif isinstance(arg, str):
        bstr = arg.encode()
    elif hasattr(arg, '__code__'):
        return hash_code(arg)
    elif isinstance(arg, functools.partial):
        return hash_partial(arg)
    elif isinstance(arg, np.ndarray):
        bstr = arg.tobytes()
    else:
        try:
            return hash(arg)
        except TypeError:
            bstr = str(arg).encode()
    m = hashlib.md5()
    m.update(bstr)
    return int.from_bytes(m.digest(), 'big')


def merge_hashes(*hashes):
    """Merge several hashes to one hash value."""
    return hash_anything(''.join([str(h) for h in hashes]))


def five_point_stencil(xdata, ydata):
    """
    Calculate the derivative dy/dx with a five point stencil.
    This algorith is only valid for equally distributed x values.

    Args:
        xdata: x values of the data points
        ydata: y values of the data points

    Returns:
        Values where the derivative was estimated and the value of the derivative at these points.

    This algorithm is only valid for values on a regular grid, for unevenly distributed
    data it is only an approximation, albeit a quite good one.

    See: https://en.wikipedia.org/wiki/Five-point_stencil
    """
    return xdata[2:-2], (
        (-ydata[4:] + 8 * ydata[3:-1] - 8 * ydata[1:-3] + ydata[:-4]) /
        (3 * (xdata[4:] - xdata[:-4]))
        )


def filon_fourier_transformation(time, correlation,
                                 frequencies=None, derivative='linear', imag=True,
                                 ):
    """
    Fourier-transformation for slow varrying functions. The filon algorithmus is
    described in detail in ref [1], ch. 3.2.3.

    Args:
        time: List of times where the correlation function was sampled.
        correlation: Values of the correlation function.
        frequencies (opt.):
            List of frequencies where the fourier transformation will be calculated.
            If None the frequencies will be choosen based on the input times.
        derivative (opt.):
            Approximation algorithmus for the derivative of the correlation function.
            Possible values are: 'linear', 'stencil' or a list of derivatives.
        imag (opt.): If imaginary part of the integral should be calculated.

    If frequencies are not explicitly given they will be evenly placed on a log scale
    in the interval [1/tmax, 0.1/tmin] where tmin and tmax are the smallest respectively
    the biggest time (greater than 0) of the provided times. The frequencies are cut off
    at high values by one decade, since the fourier transformation deviates quite strongly
    in this regime.

    References:
        [1] T. Blochowicz, Broadband dielectric spectroscopy in neat and binary
        molecular glass formers, Ph.D. thesis, Uni-versität Bayreuth (2003)
    """
    if frequencies is None:
        f_min = 1 / time[time > 0][-1]
        f_max = 0.05**(1.2 - correlation[correlation > 0][0]) / time[time > 0][0]
        frequencies = 2 * np.pi * np.logspace(
            np.log10(f_min), np.log10(f_max), num=100
        )
    frequencies.reshape(1, -1)

    if derivative is 'linear':
        derivative = (np.diff(correlation) / np.diff(time)).reshape(-1, 1)
    elif derivative is 'stencil':
        _, derivative = five_point_stencil(time, correlation)
        time = ((time[2:-1] * time[1:-2])**.5).reshape(-1, 1)
        derivative = derivative.reshape(-1, 1)
    elif np.iterable(derivative) and len(time) is len(derivative):
        derivative.reshape(-1, 1)
    else:
        raise NotImplementedError(
            'Invalid approximation method {}. Possible values are "linear", "stencil" or a list of values.'
        )
    time = time.reshape(-1, 1)

    integral = (np.cos(frequencies * time[1:]) - np.cos(frequencies * time[:-1])) / frequencies**2
    fourier = (derivative * integral).sum(axis=0)

    if imag:
        integral = 1j * (np.sin(frequencies * time[1:]) - np.sin(frequencies * time[:-1])) / frequencies**2
        fourier += (derivative * integral).sum(axis=0) + 1j * correlation[0] / frequencies

    return frequencies.reshape(-1,), fourier


def mask2indices(mask):
    """
    Return the selected indices of an array mask.
    If the mask is two-dimensional, the indices will be calculated for the second axis.

    Example:
        >>> mask2indices([True, False, True, False])
        array([0, 2])
        >>> mask2indices([[True, True, False], [True, False, True]])
        array([[0, 1], [0, 2]])
    """
    mask = np.array(mask)
    if len(mask.shape) == 1:
        indices = np.where(mask)
    else:
        indices = np.array([np.where(m) for m in mask])
    return indices


def superpose(x1, y1, x2, y2, N=100, damping=1.0):
    if x2[0] == 0:
        x2 = x2[1:]
        y2 = y2[1:]

    reg1 = x1 < x2[0]
    reg2 = x2 > x1[-1]
    x_ol = np.logspace(
        np.log10(max(x1[~reg1][0], x2[~reg2][0]) + 0.001),
        np.log10(min(x1[~reg1][-1], x2[~reg2][-1]) - 0.001),
        (sum(~reg1)+sum(~reg2))/2
    )

    def w(x):
        A = x_ol.min()
        B = x_ol.max()
        return (np.log10(B / x)/np.log10(B / A))**damping

    xdata = np.concatenate((
            x1[reg1],
            x_ol,
            x2[reg2]))
    y1_interp = interp1d(x1[~reg1], y1[~reg1])
    y2_interp = interp1d(x2[~reg2], y2[~reg2])
    ydata = np.concatenate((
        y1[x1 < x2.min()],
        w(x_ol)*y1_interp(x_ol) + (1 - w(x_ol))*y2_interp(x_ol),
        y2[x2 > x1.max()]
        ))
    return xdata, ydata


def runningmean(data, nav):
    """
    Compute the running mean of a 1-dimenional array for `nav` points.
    """
    return np.convolve(data, np.ones((nav,))/nav, mode='valid')


def coherent_sum(func, coord_a, coord_b):
    """
    Perform a coherent sum over two arrays: :math:`\frac{1}{N_a N_b}\sum_i\sum_j f(X_a[i], X_b[j])`.
    For numpy arrays ofthis is equal to::

        N, d = x.shape
        M, d = y.shape
        coherent_sum(f, x, y) == f(x.reshape(N, 1, d), x.reshape(1, M, d)).sum()

    Args:
        func: The function is called for each two items in both arrays, this should return a scalar value.
        coord_a, coord_b: The two arrays.
        normed (opt.): If the result should be normalized by the number of summands.

    """
    if isinstance(func, FunctionType):
        func = numba.jit(func, nopython=True, cache=True)
    @numba.jit(nopython=True)
    def cohsum(coord_a, coord_b):
        res = 0
        for i in range(len(coord_a)):
            for j in range(len(coord_b)):
                res += func(coord_a[i], coord_b[j])
        return res

    return cohsum(coord_a, coord_b)


@numba.jit
def norm(vec):
    return (vec**2).sum()**0.5
