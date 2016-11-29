import numpy as np
import numba
from scipy.special import legendre
from itertools import chain
from functools import reduce

from .pbc import pbc_diff, pbc_diff_numba
from .meta import annotate
from .autosave import autosave_data
from .utils import filon_fourier_transformation, norm, coherent_sum


def log_indices(first, last, num=100):
    ls = np.logspace(0, np.log10(last - first + 1), num=num)
    return np.unique(np.int_(ls) - 1 + first)


def correlation(function, frames):
    iterator = iter(frames)
    start_frame = next(iterator)
    return map(lambda f: function(start_frame, f), chain([start_frame], iterator))


def subensemble_correlation(selector_function, correlation_function=correlation):

    def c(function, frames):
        iterator = iter(frames)
        start_frame = next(iterator)
        selector = selector_function(start_frame)
        subensemble = map(lambda f: f[selector], chain([start_frame], iterator))
        return correlation_function(function, subensemble)
    return c


@autosave_data(nargs=2, kwargs_keys=(
    'index_distribution', 'correlation', 'segments', 'window', 'skip', 'average'
))
def shifted_correlation(function, frames,
                        index_distribution=log_indices, correlation=correlation,
                        segments=10, window=0.5, skip=0,
                        average=False,):
    """
    Calculate the time series for a correlation function

    The times at which the correlation is calculated are determined automatically by the
    function given as ``index_distribution``. The default is a logarithmic distribution.

    Args:
        function:   The function that should be correlated
        frames:     The coordinates of the simulation data
        index_distribution (opt.):
                    A function that returns the indices for which the timeseries
                    will be calculated
        correlation (function, opt.):
                    The correlation function
        segments (int, opt.):
                    The number of segments the time window will be shifted
        window (float, opt.):
                    The fraction of the simulation the time series will cover
        skip (float, opt.):
                    The fraction of the trajectory that will be skipped at the beginning.
    Returns:
        tuple:
            A list of length N that contains the indices of the frames at which
            the time series was calculated and a numpy array of shape (segments, N)
            that holds the (non-avaraged) correlation data

    Example:
        Calculating the mean square displacement of a coordinates object named ``coords``:

        >>> indices, data = shifted_correlation(msd, coords)
    """
    assert window + skip < 1
    start_frames = np.int_(np.linspace(len(frames)*skip, len(frames) * (1 - window - skip), num=segments, endpoint=False))
    num_frames = int(len(frames) * (window))

    idx = index_distribution(0, num_frames)

    def correlate(start_frame):
        shifted_idx = idx + start_frame
        return correlation(function, map(frames.__getitem__, shifted_idx))

    result = np.array([list(correlate(start_frame)) for start_frame in start_frames])
    if average:
        result = result.mean(axis=0)
    times = np.array([frames[i].time for i in idx]) - frames[0].time
    return times, result


def msd(start, frame, box=None):
    """
    Mean square displacement
    """
    vec = pbc_diff(start, frame, box)
    return (vec ** 2).sum(axis=1).mean()


def isf(start, frame, q, box=None):
    """
    Incoherent intermediate scattering function. To specify q, use
    water_isf = functools.partial(isf, q=22.77) # q has the value 22.77 nm^-1

    :param q: length of scattering vector
    """
    vec = pbc_diff(start, frame, box)  # start-frame
    distance = (vec ** 2).sum(axis=1) ** .5
    return np.sinc(distance * q / np.pi).mean()


def rotational_autocorrelation(onset, frame, order=2):
    """
    Compute the rotaional autocorrelation of the legendre polynamial for the given vectors.

    Args:
        onset, frame: CoordinateFrames of vectors
        order (opt.): Order of the legendre polynomial.

    Returns:
        Skalar value of the correltaion function.
    """
    scalar_prod = (onset * frame).sum(axis=-1)
    poly = legendre(order)
    return poly(scalar_prod).mean()


@annotate.untested
def oaf(start, frame):
    """
    Orientation autocorrelation function. start and frame must be connection vectors, not absolute coordinates. Use for
    example oaf_indexed to define connection vectors.

    :param start:
    :param frame:
    :return:
    """
    vec_start_norm = np.norm(start)
    vec_frame_norm = np.norm(frame)

    dot_prod = (start * frame).sum(axis=1) / (vec_start_norm, vec_frame_norm)
    return (3 * dot_prod**2 - 1).mean() / 2.0


def oaf_indexed(index_from, index_to):
    """
    Returns a OAF correlation function. Example
    oaf_indexed(t[:,1] == 'C', t[:,1] == 'O')
    :param index_from:
    :param index_to:
    :return:
    """
    return lambda start, frame: oaf(start[index_to] - start[index_from],
                                    frame[index_to] - frame[index_from])


@annotate.unfinished
def van_hove(start, end, bins, box=None):
    vec = pbc_diff(start, end, box)
    delta_r = ((vec)**2).sum(axis=1) ** .5

    return 1 / len(start) * np.histogram(delta_r, bins)[0]


def overlap(onset, frame, crds_tree, radius):
    """
    Compute the overlap with a reference configuration defined in a CoordinatesTree.

    Args:
        onset: Initial frame, this is only used to get the frame index
        frame: The current configuration
        crds_tree: A CoordinatesTree of the reference configurations
        radius: The cutoff radius for the overlap

    This function is intended to be used with :func:`shifted_correlation`.
    As usual the first two arguments are used internally and the remaining ones
    should be defined with :func:`functools.partial`.

    If the overlap of a subset of the system should be calculated, this has to be
    defined through a selection of the reference configurations in the CoordinatesTree.

    Example:
        >>> shifted_correlation(
        ...     partial(overlap, crds_tree=CoordinatesTree(traj), radius=0.11),
        ...     traj
        ... )
    """
    tree = crds_tree[onset.step]
    return (tree.query(frame)[0] <= radius).sum() / tree.n


def susceptibility(time, correlation, **kwargs):
    """
    Calculate the susceptibility of a correlation function.

    Args:
        time: Timesteps of the correlation data
        correlation: Value of the correlation function
        **kwargs (opt.):
            Additional keyword arguments will be passed to :func:`filon_fourier_transformation`.
    """
    frequencies, fourier = filon_fourier_transformation(time, correlation, imag=False, **kwargs)
    return frequencies, frequencies * fourier


def coherent_scattering_function(onset, frame, q):
    """
    Calculate the coherent scattering function.
    """
    box = onset.box.diagonal()
    dimension = len(box)

    @numba.jit(nopython=True, cache=True)
    def scfunc(x, y):
        sqdist = 0
        for i in range(dimension):
            d = x[i] - y[i]
            if d > box[i] / 2:
                d -= box[i]
            if d < -box[i] / 2:
                d += box[i]
            sqdist += d**2
        x = sqdist**0.5 * q
        return np.sin(x) / x

    return coherent_sum(scfunc, onset.pbc, frame.pbc) / len(onset)
