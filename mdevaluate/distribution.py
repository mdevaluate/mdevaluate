import numpy as np
from .coordinates import pbc_diff
from .atoms import next_neighbors


def rotate_axis(coords, axis):
    """
    Rotate a set of coordinates to a given axis.
    """
    axis = np.array(axis) / np.linalg.norm(axis)
    zaxis = np.array([0, 0, 1])
    if (axis == zaxis).sum() == 3:
        return coords
    rotation_axis = np.cross(axis, zaxis)
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)

    theta = np.arccos(axis @ zaxis / np.linalg.norm(axis))

    # return theta/pi, rotation_axis

    ux, uy, uz = rotation_axis
    cross_matrix = array([
        [0, -uz, uy],
        [uz,  0, -ux],
        [-uy, ux, 0]
    ])
    rotation_matrix = np.cos(theta) * np.identity(len(axis)) \
        + (1 - np.cos(theta)) * rotation_axis.reshape(-1, 1) @ rotation_axis.reshape(1, -1) \
        + np.sin(theta) * cross_matrix

    if len(coords.shape) == 2:
        rotated = np.array([rotation_matrix @ xyz for xyz in coords])
    else:
        rotated = rotation_matrix @ coords
    return rotated


def polar_coordinates(x, y):
    """Convert cartesian to polar coordinates."""
    radius = np.sqrt(x**2 + y**2)
    phi = np.arctan(y / x)
    return radius, phi


def spherical_coordinates(x, y, z):
    """Convert cartesian to spherical coordinates."""
    radius, phi = polar_coordinates(x,y)
    theta = np.arccos(z / radius)
    return radius, phi, theta


def time_average(function, coordinates, pool=None, verbose=False):
    """
    Compute the time average of a function.

    Args:
        function:
            The function that will be averaged, it has to accept exactly one argument
            which is the current atom set.
        coordinates: The coordinates object of the simulation
        pool (multiprocessing.Pool, opt.):
            A multiprocessing pool which will be used for cocurrent calculation of the
            averaged function.
        verbose (bool, opt.): If

    """
    coordinate_iter = iter(coordinates)
    number_of_averages = 0
    result = 0

    if pool is not None:
        _map = pool.imap
    else:
        _map = map

    evaluated = _map(function, coordinate_iter)

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
    volume = 2 * np.pi * (bins[1:]**2 - bins[:-1]**2) * height
    return hist / volume

def shell_density(atoms, shell_radius, bins, pbc_box=None, shell_thickness=0.5, symmetry_axis=(0,0,1), origin=(0,0,0)):
    """
    Compute the density distribution on a cylindrical shell.

    Args:
        atoms: The coordinates of the atoms
        shell_radius: Inner radius of the shell
        bins: Histogramm bins, this has to be a two-dimensional list of bins: [angle, z]
        pbc_box (opt.): Tuple of box lengths in x,y,z direction which is used to set atoms
            back to the pbc-box
        shell_thickness (opt.): Thicknes of the shell, default is 0.5
        symmetry_axis (opt.): The symmtery axis of the pore, the coordinates will be
            rotated such that this axis is the z-axis
        origin (opt.): Origin of the pore, the coordinates will be moved such that this
            is the new origin.

    Returns:
        Two-dimensional density distribution of the atoms in the defined shell.
    """
    if pbc_z is not None:
        atoms[:] %= pbc_box
    cartesian = rotate_axis(atoms-origin, symmetry_axis)
    radius, theta = polar_coordinates(cartesian[:,0], cartesian[:,1])
    shell_indices = (shell_radius <= radius)&(radius <= shell_radius + shell_thickness)
    hist = numpy.histogram2d(theta[shell_indices], cartesian[shell_indices, 2], bins)[0]

    return hist
