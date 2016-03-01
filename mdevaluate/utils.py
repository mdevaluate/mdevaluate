"""
Collection of utility functions.
"""
import numpy as np


def hash_anything(arg):
    """Return a hash value for the current state of any argument."""
    try:
        return hash(arg)
    except TypeError:
        return hash(str(arg))


def merge_hashes(*hashes):
    """Merge several hashes to one hash value."""
    return hash(''.join([str(h) for h in hashes]))


def five_point_stencil(xdata, ydata):
    """
    Calculate the derivative dy/dx with a five point stencil.
    This algorith is only valid for equally distributed x values.

    Args:
        xdata: x values of the data points
        ydata: y values of the data points

    Returns:
        Values where the derivative was estimated and the value of the derivative at these points.

    See: https://en.wikipedia.org/wiki/Five-point_stencil
    """
    return xdata[1:-1], (
        (-ydata[3:] + 8 * ydata[2:-1] - 8 * ydata[1:-2] + ydata[:-3]) /
        (12 * (xdata[2:-1] - xdata[1:-2]))
        )


def filon_fourier_transformation(time, correlation, frequencies=None, interpolation='linear'):
    """
    Fourier-transformation for slow varrying functions. The filon algorithmus is
    described in detail in [1].

    Args:
        time: List of times where the correlation function was sampled.
        correlatio: Values of the correlation function.
        frequencies (opt.):
            List of frequencies where the fourier transformation will be calculated.
            If None the frequencies will be choosen based on the input times.
        interpolation (opt.):
            Interpolation algorithmus for the derivative of the correlation function.
            Possible values are: 'linear' or 'stencil'.

    Reference:
        [1] T. Blochowicz, Broadband dielectric spectroscopy in neat and binary
        molecular glass formers, Ph.D. thesis, Uni-versität Bayreuth (2003)
    """
    if frequencies is None:
        frequencies = 2*np.pi*np.logspace(
            np.log10(1 / time[-1]), np.log10(1 / time[0]), num=100
        )
    frequencies.reshape(1, -1)

    if interpolation is 'linear':
        derivative = (np.diff(correlation) / np.diff(time)).reshape(-1, 1)
    elif interpolation is 'stencil':
        time, derivative = five_point_stencil(time, correlation).reshape(-1, 1)
    else:
        raise NotImplementedError(
            'Invalid interpolation method {}. Possible values are "linear" or "stencil" '
            )
    time = time.reshape(-1, 1)

    integral = (np.cos(frequencies * time[1:]) - np.cos(frequencies * time[:-1])) / frequencies**2
    fourier = (derivative * integral).sum(axis=0) / derivative.size

    return frequencies.reshape(-1,), fourier
