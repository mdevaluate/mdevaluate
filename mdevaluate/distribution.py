import warnings

import numpy as np

from .coordinates import pbc_diff, rotate_axis, polar_coordinates, spherical_coordinates
from .atoms import next_neighbors
from .autosave import autosave_data


@autosave_data(nargs=2, kwargs_keys=('coordinates_b',))
def time_average(function, coordinates, coordinates_b=None, pool=None, verbose=False):
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
        verbose (bool, opt.): Be verbose about the progress

    """
    if pool is not None:
        _map = pool.imap
    else:
        _map = map

    number_of_averages = 0
    result = 0

    if coordinates_b is not None:
        coordinate_iter = (iter(coordinates), iter(coordinates_b))
    else:
        coordinate_iter = (iter(coordinates),)

    evaluated = _map(function, *coordinate_iter)

    for ev in evaluated:
        number_of_averages += 1
        result += ev
        if verbose and number_of_averages % 100 == 0:
            print(number_of_averages)

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


def radial_pair_distribution(atoms, bins, box=None):
    number_of_atoms = len(atoms)
    # Calculate the upper triangle of differences
    indices = np.triu_indices(number_of_atoms, k=1)

    # vector between atoms
    diff = pbc_diff(atoms[indices[0]], atoms[indices[1]], box)

    dist = (diff**2).sum(axis=1)**.5
    volume = 4 / 3 * np.pi * (bins[1:]**3 - bins[:-1]**3)
    hist = np.histogram(dist, bins)[0]
    density = (number_of_atoms - 1) / np.prod(box)

    return hist / volume / density * (2 / (number_of_atoms - 1))


def rdf(atoms_a, atoms_b=None, bins=None, box=None):
    """
    Compute the radial pair distribution of one or two sets of atoms.

    .. math::
       g_{AB}(r) = \\frac{1}{\langle \\rho_B\\rangle N_A}\sum\limits_{i\in A}^{N_A}
       \sum\limits_{j\in B}^{N_B}\\frac{\delta(r_{ij} -r)}{4\pi r^2}

    For use with :func:`time_average`, define bins through the use of :func:`~functools.partial`,
    the atom sets are passed to :func:`time_average`, if a second set of atoms should be used
    specify it as ``coordinates_b`` and it will be passed to this function.

    Args:
        atoms_a: First set of atoms, used internally
        atoms_b (opt.): Second set of atoms, used internally
        bins: Bins of the radial distribution function
        box (opt.): Simulations box, if not specified this is taken from ``atoms_a.box``
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
        if nr_a > nr_b:
            # Make sure atoms_a is always the smaller set.
            atoms_a, atoms_b = atoms_b, atoms_a
        fill = np.empty((abs(nr_a - nr_b), dim))
        fill.fill(np.nan)
        atoms_a = np.concatenate((atoms_a, fill))
        nr_of_atoms = max(nr_a, nr_b)
        indices = np.triu_indices(nr_of_atoms, k=0)

    with warnings.catch_warnings(record=True):
        diff = pbc_diff(atoms_a[indices[0]], atoms_b[indices[1]], box)
    diff = diff[~np.isnan(diff).any(axis=1)]
    dist = (diff**2).sum(axis=1)**0.5
    volume = 4 / 3 * np.pi * (bins[1:]**3 - bins[:-1]**3)
    hist, _ = np.histogram(dist, bins)
    density = len(diff) / np.prod(box)

    return hist / volume / density


def distance_distribution(atoms, bins):
    connection_vectors = atoms[:-1, :] - atoms[1:, :]
    connection_lengths = (connection_vectors**2).sum(axis=1)**.5
    return np.histogram(connection_lengths, bins)[0]


def tetrahedral_order(atoms):
    indices = next_neighbors(atoms, number_of_neighbors=4)
    neighbors = atoms[indices]
    neighbors_1, neighbors_2, neighbors_3, neighbors_4 = \
        neighbors[:, 0, :], neighbors[:, 1, :], neighbors[:, 2, :], neighbors[:, 3, :]

    # Connection vectors
    neighbors_1 -= atoms
    neighbors_2 -= atoms
    neighbors_3 -= atoms
    neighbors_4 -= atoms

    # Normed Connection vectors
    neighbors_1 /= ((neighbors_1**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0], 1)
    neighbors_2 /= ((neighbors_2**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0], 1)
    neighbors_3 /= ((neighbors_3**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0], 1)
    neighbors_4 /= ((neighbors_4**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0], 1)

    a_1_2 = ((neighbors_1 * neighbors_2).sum(axis=1) + 1 / 3)**2
    a_1_3 = ((neighbors_1 * neighbors_3).sum(axis=1) + 1 / 3)**2
    a_1_4 = ((neighbors_1 * neighbors_4).sum(axis=1) + 1 / 3)**2

    a_2_3 = ((neighbors_2 * neighbors_3).sum(axis=1) + 1 / 3)**2
    a_2_4 = ((neighbors_2 * neighbors_4).sum(axis=1) + 1 / 3)**2

    a_3_4 = ((neighbors_3 * neighbors_4).sum(axis=1) + 1 / 3)**2

    q = 1 - 3 / 8 * (a_1_2 + a_1_3 + a_1_4 + a_2_3 + a_2_4 + a_3_4)

    return q


def tetrahedral_order_distribution(atoms, bins):
    Q = tetrahedral_order(atoms)
    return np.histogram(Q, bins=bins)[0]


def radial_density(atoms, bins, symmetry_axis=(0, 0, 1), origin=(0, 0, 0), height=1):
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
    """
    cartesian = rotate_axis(atoms - origin, symmetry_axis)
    radius, _ = polar_coordinates(cartesian[:, 0], cartesian[:, 1])
    hist = np.histogram(radius, bins=bins)[0]
    volume = np.pi * (bins[1:]**2 - bins[:-1]**2) * height
    return hist / volume


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
