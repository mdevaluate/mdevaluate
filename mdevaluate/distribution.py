import warnings
import numpy as np

from .coordinates import rotate_axis, polar_coordinates, spherical_coordinates
from .atoms import next_neighbors
from .autosave import autosave_data
from .meta.annotate import deprecated
from .utils import runningmean
from .pbc import pbc_diff, pbc_points
from .logging import logger
from scipy import spatial


@autosave_data(nargs=2, kwargs_keys=('coordinates_b',), version='time_average-1')
def time_average(function, coordinates, coordinates_b=None, pool=None):
    """
    Compute the time average of a function.

    Args:
        function:
            The function that will be averaged, it has to accept exactly one argument
            which is the current atom set
        coordinates: The coordinates object of the simulation
        pool (multiprocessing.Pool, opt.):
            A multiprocessing pool which will be used for cocurrent calculation of the
            averaged function

    """
    if pool is not None:
        _map = pool.imap
    else:
        _map = map

    number_of_averages = 0
    result = 0

    if coordinates_b is not None:
        if coordinates._slice != coordinates_b._slice:
            logger.warning("Different slice fro coordinates and coordinates_b.")
        coordinate_iter = (iter(coordinates), iter(coordinates_b))
    else:
        coordinate_iter = (iter(coordinates),)

    evaluated = _map(function, *coordinate_iter)

    for ev in evaluated:
        number_of_averages += 1
        result += ev
        if number_of_averages % 100 == 0:
            logger.debug('time_average: %d', number_of_averages)

    return result / number_of_averages


def time_histogram(function, coordinates, bins, hist_range, pool=None):
    coordinate_iter = iter(coordinates)

    if pool is not None:
        _map = pool.imap
    else:
        _map = map

    evaluated = _map(function, coordinate_iter)

    results = []
    hist_results = []
    for num, ev in enumerate(evaluated):
        results.append(ev)

        if num % 100 == 0 and num > 0:
            print(num)
        r = np.array(results).T
        for i, row in enumerate(r):
            histo, _ = np.histogram(row, bins=bins, range=hist_range)
            if len(hist_results) <= i:
                hist_results.append(histo)
            else:
                hist_results[i] += histo
            results = []
    return hist_results


def rdf(atoms_a, atoms_b=None, bins=None, box=None, chunksize=50000, returnx=False, **kwargs):
    r"""
    Compute the radial pair distribution of one or two sets of atoms.

    .. math::
       g_{AB}(r) = \frac{1}{\langle \rho_B\rangle N_A}\sum\limits_{i\in A}^{N_A}
       \sum\limits_{j\in B}^{N_B}\frac{\delta(r_{ij} -r)}{4\pi r^2}

    For use with :func:`time_average`, define bins through the use of :func:`~functools.partial`,
    the atom sets are passed to :func:`time_average`, if a second set of atoms should be used
    specify it as ``coordinates_b`` and it will be passed to this function.

    Args:
        atoms_a: First set of atoms, used internally
        atoms_b (opt.): Second set of atoms, used internally
        bins: Bins of the radial distribution function
        box (opt.): Simulations box, if not specified this is taken from ``atoms_a.box``
        chunksize (opt.):
            For large systems (N > 1000) the distaces have to be computed in chunks so the arrays
            fit into memory, this parameter controlls the size of these chunks. It should be
            as large as possible, depending on the available memory.
        returnx (opt.): If True the x ordinate of the histogram is returned.
    """
    assert bins is not None, 'Bins of the pair distribution have to be defined.'
    if box is None:
        box = atoms_a.box.diagonal()
    if atoms_b is None:
        atoms_b = atoms_a
        nr_of_atoms = len(atoms_a)
        indices = np.triu_indices(nr_of_atoms, k=1)
    else:
        nr_a, dim = atoms_a.shape
        nr_b, dim = atoms_b.shape
        indices = np.array([(i, j) for i in range(nr_a) for j in range(nr_b)]).T

    # compute the histogram in chunks for large systems
    hist = 0
    nr_of_samples = 0
    for chunk in range(0, len(indices[0]), chunksize):
        sl = slice(chunk, chunk + chunksize)
        diff = pbc_diff(atoms_a[indices[0][sl]], atoms_b[indices[1][sl]], box)
        dist = (diff**2).sum(axis=1)**0.5
        nr_of_samples += len(dist)
        hist += np.histogram(dist, bins)[0]

    volume = 4 / 3 * np.pi * (bins[1:]**3 - bins[:-1]**3)
    density = nr_of_samples / np.prod(box)
    res = hist / volume / density
    if returnx:
        return np.vstack((runningmean(bins, 2), res))
    else:
        return res


def pbc_tree_rdf(atoms_a, atoms_b=None, bins=None, box=None, exclude=0, returnx=False, **kwargs):
    if box is None:
        box = atoms_a.box.diagonal()
    all_coords = pbc_points(pbc_diff(atoms_b,box=box), box, thickness=np.amax(bins)+0.1, center=0)
    to_tree = spatial.cKDTree(all_coords)
    dist = to_tree.query(pbc_diff(atoms_a,box=box),k=len(atoms_b), distance_upper_bound=np.amax(bins)+0.1)[0].flatten()
    dist = dist[dist < np.inf]
    hist = np.histogram(dist, bins)[0]
    volume = 4/3*np.pi*(bins[1:]**3-bins[:-1]**3)
    res = (hist) * np.prod(box) / volume / len(atoms_a) / (len(atoms_b)-exclude)
    if returnx:
        return np.vstack((runningmean(bins, 2), res))
    else:
        return res

def pbc_spm_rdf(atoms_a, atoms_b=None, bins=None, box=None, exclude=0, returnx=False, **kwargs):
    if box is None:
        box = atoms_a.box.diagonal()
    all_coords = pbc_points(pbc_diff(atoms_b,box=box), box, thickness=np.amax(bins)+0.1, center=0)
    to_tree = spatial.cKDTree(all_coords)
    if all_coords.nbytes/1024**3 * len(atoms_a) < 2:
        from_tree = spatial.cKDTree(pbc_diff(atoms_a,box=box))
        dist = to_tree.sparse_distance_matrix(from_tree, max_distance=np.amax(bins)+0.1, output_type='ndarray')
        dist = np.asarray(dist.tolist())[:,2]
        hist = np.histogram(dist, bins)[0]
    else:
        chunksize = int(2 * len(atoms_a) / (all_coords.nbytes/1024**3 * len(atoms_a)))
        hist = 0
        for chunk in range(0, len(atoms_a), chunksize):
            sl = slice(chunk, chunk + chunksize)
            from_tree = spatial.cKDTree(pbc_diff(atoms_a[sl],box=box))
            dist = to_tree.sparse_distance_matrix(from_tree, max_distance=np.amax(bins)+0.1, output_type='ndarray')
            dist = np.asarray(dist.tolist())[:,2]
            hist += np.histogram(dist, bins)[0]

    volume = 4/3*np.pi*(bins[1:]**3-bins[:-1]**3)
    res = (hist) * np.prod(box) / volume / len(atoms_a) / (len(atoms_b)-exclude)
    if returnx:
        return np.vstack((runningmean(bins, 2), res))
    else:
        return res

@autosave_data(nargs=2, kwargs_keys=('to_coords','times'))
def fast_averaged_rdf(from_coords, bins, to_coords=None, times=10, exclude=0, **kwargs):
    if to_coords is None:
        to_coords = from_coords
        exclude = 1
    # first find timings for the different rdf functions
    import time
    # only consider sparse matrix for this condition
    if (len(from_coords[0])*len(to_coords[0]) <= 3000 * 2000 ) & (len(from_coords[0])/len(to_coords[0]) > 5 ):
        funcs = [rdf, pbc_tree_rdf, pbc_spm_rdf]
    else:
        funcs = [rdf, pbc_tree_rdf]
    timings = []
    for f in funcs:
        start = time.time()
        f(from_coords[0], atoms_b=to_coords[0], bins=bins, box=np.diag(from_coords[0].box))
        end = time.time()
        timings.append(end-start)
    timings = np.array(timings)
    timings[0] = 2*timings[0] # statistics for the other functions is twice as good per frame
    logger.debug('rdf function timings: ' + str(timings))
    rdffunc = funcs[np.argmin(timings)]
    logger.debug('rdf function used: ' + str(rdffunc))
    if rdffunc == rdf:
        times = times*2 # duplicate times for same statistics

    frames = np.array(range(0, len(from_coords), int(len(from_coords)/times)))[:times]
    out = np.zeros(len(bins)-1)
    for j, i in enumerate(frames):
        logger.debug('multi_radial_pair_distribution: %d/%d', j, len(frames))
        out += rdffunc(from_coords[i], to_coords[i], bins, box=np.diag(from_coords[i].box), exclude=exclude)
    return out/len(frames)


def distance_distribution(atoms, bins):
    connection_vectors = atoms[:-1, :] - atoms[1:, :]
    connection_lengths = (connection_vectors**2).sum(axis=1)**.5
    return np.histogram(connection_lengths, bins)[0]


def tetrahedral_order(atoms, reference_atoms=None):
    if reference_atoms is None:
        reference_atoms = atoms
    indices = next_neighbors(reference_atoms, query_atoms=atoms, number_of_neighbors=4)
    neighbors = reference_atoms[indices]
    neighbors_1, neighbors_2, neighbors_3, neighbors_4 = \
        neighbors[:, 0, :], neighbors[:, 1, :], neighbors[:, 2, :], neighbors[:, 3, :]

    # Connection vectors
    neighbors_1 -= atoms
    neighbors_2 -= atoms
    neighbors_3 -= atoms
    neighbors_4 -= atoms

    # Normed Connection vectors
    neighbors_1 /= np.linalg.norm(neighbors_1, axis=-1).reshape(-1, 1)
    neighbors_2 /= np.linalg.norm(neighbors_2, axis=-1).reshape(-1, 1)
    neighbors_3 /= np.linalg.norm(neighbors_3, axis=-1).reshape(-1, 1)
    neighbors_4 /= np.linalg.norm(neighbors_4, axis=-1).reshape(-1, 1)

    a_1_2 = ((neighbors_1 * neighbors_2).sum(axis=1) + 1 / 3)**2
    a_1_3 = ((neighbors_1 * neighbors_3).sum(axis=1) + 1 / 3)**2
    a_1_4 = ((neighbors_1 * neighbors_4).sum(axis=1) + 1 / 3)**2

    a_2_3 = ((neighbors_2 * neighbors_3).sum(axis=1) + 1 / 3)**2
    a_2_4 = ((neighbors_2 * neighbors_4).sum(axis=1) + 1 / 3)**2

    a_3_4 = ((neighbors_3 * neighbors_4).sum(axis=1) + 1 / 3)**2

    q = 1 - 3 / 8 * (a_1_2 + a_1_3 + a_1_4 + a_2_3 + a_2_4 + a_3_4)

    return q


def tetrahedral_order_distribution(atoms, reference_atoms=None, bins=None):
    assert bins is not None, 'Bin edges of the distribution have to be specified.'
    Q = tetrahedral_order(atoms, reference_atoms=reference_atoms)
    return np.histogram(Q, bins=bins)[0]


def radial_density(atoms, bins, symmetry_axis=(0, 0, 1), origin=(0, 0, 0), height=1, returnx=False):
    """
    Calculate the radial density distribution.

    This function is meant to be used with time_average.

    Args:
        atoms:
            Set of coordinates.
        bins:
            Bin specification that is passed to numpy.histogram. This needs to be
            a list of bin edges if the function is used within time_average.
        symmetry_axis (opt.):
            Vector of the symmetry axis, around which the radial density is calculated,
            default is z-axis.
        origin (opt.):
            Origin of the rotational symmetry, e.g. center of the pore.
        height (opt.):
            Height of the pore, necessary for correct normalization of the density.
        returnx (opt.):
            If True, the x ordinate of the distribution is returned.
    """
    cartesian = rotate_axis(atoms - origin, symmetry_axis)
    radius, _ = polar_coordinates(cartesian[:, 0], cartesian[:, 1])
    hist = np.histogram(radius, bins=bins)[0]
    volume = np.pi * (bins[1:]**2 - bins[:-1]**2) * height
    res = hist / volume
    if returnx:
        return np.vstack((runningmean(bins, 2), res))
    else:
        return res


def shell_density(atoms, shell_radius, bins, shell_thickness=0.5,
                  symmetry_axis=(0, 0, 1), origin=(0, 0, 0)):
    """
    Compute the density distribution on a cylindrical shell.

    Args:
        atoms: The coordinates of the atoms
        shell_radius: Inner radius of the shell
        bins: Histogramm bins, this has to be a two-dimensional list of bins: [angle, z]
        shell_thickness (opt.): Thicknes of the shell, default is 0.5
        symmetry_axis (opt.): The symmtery axis of the pore, the coordinates will be
            rotated such that this axis is the z-axis
        origin (opt.): Origin of the pore, the coordinates will be moved such that this
            is the new origin.

    Returns:
        Two-dimensional density distribution of the atoms in the defined shell.
    """
    cartesian = rotate_axis(atoms-origin, symmetry_axis)
    radius, theta = polar_coordinates(cartesian[:, 0], cartesian[:, 1])
    shell_indices = (shell_radius <= radius) & (radius <= shell_radius + shell_thickness)
    hist = np.histogram2d(theta[shell_indices], cartesian[shell_indices, 2], bins)[0]

    return hist


def spatial_density(atoms, bins, weights=None):
    """
    Compute the spatial density distribution.
    """
    density, _ = np.histogramdd(atoms, bins=bins, weights=weights)
    return density


def mixing_ratio_distribution(atoms_a, atoms_b, bins_ratio, bins_density,
                              weights_a=None, weights_b=None, weights_ratio=None):
    """
    Compute the distribution of the mixing ratio of two sets of atoms.
    """

    density_a, _ = time_average
    density_b, _ = np.histogramdd(atoms_b, bins=bins_density, weights=weights_b)
    mixing_ratio = density_a/(density_a + density_b)
    good_inds = (density_a != 0) & (density_b != 0)
    hist, _ = np.histogram(mixing_ratio[good_inds], bins=bins_ratio, weights=weights_ratio)
    return hist


def next_neighbor_distribution(atoms, reference=None, number_of_neighbors=4, bins=None, normed=True):
    """
    Compute the distribution of next neighbors with the same residue name.
    """
    assert bins is not None, 'Bins have to be specified.'
    if reference is None:
        reference = atoms
    nn = next_neighbors(reference, query_atoms=atoms, number_of_neighbors=number_of_neighbors)
    resname_nn = reference.residue_names[nn]
    count_nn = (resname_nn == atoms.residue_names.reshape(-1, 1)).sum(axis=1)
    return np.histogram(count_nn, bins=bins, normed=normed)[0]
