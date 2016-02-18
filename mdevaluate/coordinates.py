import numpy as np
from functools import wraps
from copy import copy
from .atoms import AtomSubset
from .pbc import pbc_diff
from .utils import hash_anything as _hash, merge_hashes


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
    phi = np.arctan2(y, x)
    return radius, phi


def spherical_coordinates(x, y, z):
    """Convert cartesian to spherical coordinates."""
    radius, phi = polar_coordinates(x,y)
    theta = np.arccos(z / radius)
    return radius, phi, theta


class Coordinates:
    """
    Coordinates represent trajectory data, which is used for evaluation functions.

    Atoms may be selected by specifing a atom_subset or a atom_filter.
    """

    def __init__(self, frames, atom_filter=None, atom_subset: AtomSubset=None):
        self.frames = frames
        self._slice = slice(0, len(self.frames))
        assert atom_filter is None or atom_subset is None, "Cannot use both: subset and filter"

        if atom_filter is not None:
            self.atom_filter = atom_filter
        elif atom_subset is not None:
            self.atom_filter = atom_subset.selection
            self.atom_subset = atom_subset
            self.atoms = atom_subset.atoms
        else:
            self.atom_filter = np.ones(shape=(len(frames[0].coordinates),), dtype=bool)

    def slice(self, slice):
        for i in range(len(self))[slice]:
            yield self[i]

    def __iter__(self):
        return self.slice(self._slice)

    def __getitem__(self, item):
        if isinstance(item, slice):
            sliced = copy(self)
            sliced._slice = item
            return sliced
        else:
            try:
                if self.atom_filter is not None:
                    return self.frames.__getitem__(item).coordinates[self.atom_filter]
                else:
                    return self.frames.__getitem__(item).coordinates
            except EOFError:
                raise IndexError

    def __len__(self):
        return len(self.frames)

    def __hash__(self):
        return merge_hashes(_hash(self.frames), _hash(self.atom_filter), _hash(str(self._slice)))


class MeanCoordinates(Coordinates):

    def __init__(self, frames, atom_filter=None, mean=1):
        super().__init__(frames, atom_filter)
        self.mean = mean
        assert mean >= 1, "Mean must be positive"

    def __getitem__(self, item):
        frame = super().__getitem__(item)
        for i in range(item + 1, item + self.mean):
            frame += super().__getitem__(i)

        return frame / self.mean

    def len(self):
        return len(super() - self.mean + 1)


class CoordinatesMap:

    def __init__(self, coordinates, function):
        self.coordinates = coordinates
        self.function = function

    def __iter__(self):
        for frame in self.coordinates:
            yield self.function(frame)

    def __getitem__(self, item):
        if isinstance(item, slice):
            sliced = copy(self)
            sliced.coordinates = self.coordinates[item]
            return sliced
        else:
            return self.function(self.coordinates.__getitem__(item))

    def __len__(self):
        return len(self.coordinates.frames)

    def __hash__(self):
        return merge_hashes(_hash(self.coordinates), _hash(self.function.__code__))

def map_coordinates(func):
    @wraps(func)
    def wrapped(coordinates, *args, **kwargs):
        return CoordinatesMap(coordinates, lambda x: func(x, *args, **kwargs))
    return wrapped


@map_coordinates
def centers_of_mass(c, masses):
    """

    A- 1
    B- 2
    A- 1
    C  3
    A-
    B-
    A-
    C
    A-
    B-
    A-
    C


    Example:
    rd = XTCReader('t.xtc')
    coordinates = Coordinates(rd)
    com = centers_of_mass(coordinates, (1.0, 2.0, 1.0, 3.0))

    """
    # At first, regroup our array
    number_of_masses = len(masses)
    number_of_coordinates, number_of_dimensions = c.shape
    number_of_new_coordinates = number_of_coordinates // number_of_masses
    grouped_masses = c.reshape(number_of_new_coordinates, number_of_masses, number_of_dimensions)

    return np.average(grouped_masses, axis=1, weights=masses)


@map_coordinates
def pore_coordinates(coordinates, origin, sym_axis='z'):
    """
    Map coordinates of a pore simulation so the pore has cylindrical symmetry.

    Args:
        coordinates: Coordinates of the simulation
        origin: Origiin of the pore which will be the coordinates origin after mapping
        sym_axis (opt.): Symmtery axis of the pore, may be a literal direction
            'x', 'y' or 'z' or an array of shape (3,)
    """
    if sym_axis in ('x', 'y', 'z'):
        rot_axis = np.zeros(shape=(3,))
        rot_axis[['x','y','z'].index(sym_axis)] = 1
    else:
        rot_axis = sym_axis

    return rotate_axis(coordinates - origin, rot_axis)
