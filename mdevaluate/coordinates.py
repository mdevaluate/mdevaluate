from functools import wraps
from functools import partial, lru_cache
from copy import copy

import numpy as np
from scipy.spatial import cKDTree, KDTree

from .atoms import AtomSubset
from .pbc import pbc_diff, whole
from .utils import hash_anything as _hash, merge_hashes, mask2indices


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
    cross_matrix = np.array([
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


def spherical_radius(frame, origin=None):
    """
    Transform a frame of cartesian coordinates into the sperical radius.
    If origin=None the center of the box is taken as the coordinates origin.
    """
    if origin is None:
        origin = frame.box.diagonal() / 2
    return ((frame - origin)**2).sum(axis=-1)**0.5


def polar_coordinates(x, y):
    """Convert cartesian to polar coordinates."""
    radius = (x**2 + y**2)**0.5
    phi = np.arctan2(y, x)
    return radius, phi


def spherical_coordinates(x, y, z):
    """Convert cartesian to spherical coordinates."""
    xy, phi = polar_coordinates(x, y)
    theta = np.arccos(z / xy)
    radius = (x**2 + y**2 + z**2)**0.5
    return radius, phi, theta


def radial_selector(frame, coordinates, rmin, rmax):
    """
    Return a selection of all atoms with radius in the interval [rmin, rmax].
    """
    crd = coordinates[frame.step]
    rad, _ = polar_coordinates(crd[:, 0], crd[:, 1])
    selector = (rad >= rmin) & (rad <= rmax)
    return mask2indices(selector)


def spatial_selector(frame, transform, rmin, rmax):
    """
    Select a subset of atoms which have a radius between rmin and rmax.
    Coordinates are filtered by the condition::

      rmin <= transform(frame) <= rmax

    Args:
        frame: The coordinates of the actual trajectory
        transform:
            A function that transforms the coordinates of the frames into
            the one-dimensional spatial coordinate (e.g. radius).
        rmin: Minimum value of the radius
        rmax: Maximum value of the radius
    """
    r = transform(frame)
    selector = (rmin <= r) & (rmax >= r)
    return mask2indices(selector)


class CoordinateFrame(np.ndarray):
    @property
    def box(self):
        return np.array(self.coordinates.frames[self.step].box)

    @property
    def time(self):
        return self.coordinates.frames[self.step].time

    @property
    def masses(self):
        return self.coordinates.atoms.masses[self.coordinates.atom_subset.selection]

    @property
    def charges(self):
        return self.coordinates.atoms.charges[self.coordinates.atom_subset.selection]

    @property
    def residue_ids(self):
        return self.coordinates.atom_subset.residue_ids

    @property
    def residue_names(self):
        return self.coordinates.atom_subset.residue_names

    @property
    def atom_names(self):
        return self.coordinates.atom_subset.atom_names

    @property
    def whole(self):
        return whole(self)

    def __new__(subtype, shape, dtype=float, buffer=None, offset=0, strides=None, order=None,
                coordinates=None, step=None, box=None):
        obj = np.ndarray.__new__(subtype, shape, dtype, buffer, offset, strides)

        obj.coordinates = coordinates
        obj.step = step
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self.coordinates = getattr(obj, 'coordinates', None)
        self.step = getattr(obj, 'step', None)


class Coordinates:
    """
    Coordinates represent trajectory data, which is used for evaluation functions.

    Atoms may be selected by specifing a atom_subset or a atom_filter.
    """

    def __init__(self, frames, atom_filter=None, atom_subset: AtomSubset=None, caching=False):
        """
        Args:
            frames: The trajectory reader
            atom_filter (opt.): A mask which selects a subset of the system
            atom_subset (opt.): A AtomSubset that selects a subset of the system
            caching (opt.):
                If frames should be cached. If no bool is given, the value will be used
                as the maxsize of lru_cache, which can be a number or None. Use None to
                never discard any frame.

        Note:
            The caching in Coordinates is deprecated, use the CachedReader or the function open
            from the reader module instead.
        """
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

        if caching:
            if isinstance(caching, bool):
                self.get_frame = lru_cache(maxsize=128)(self.get_frame)
            else:
                self.get_frame = lru_cache(maxsize=caching)(self.get_frame)

    def slice(self, slice):
        for i in range(len(self))[slice]:
            yield self[i]

    def get_frame(self, fnr):
        """Returns the fnr-th frame."""
        try:
            if self.atom_filter is not None:
                frame = self.frames[fnr].coordinates[self.atom_filter].view(CoordinateFrame)
            else:
                frame = self.frames.__getitem__(fnr).coordinates.view(CoordinateFrame)
            frame.coordinates = self
            frame.step = fnr
            return frame
        except EOFError:
            raise IndexError

    def clear_cache(self):
        """Clears the frame cache, if it is enabled."""
        if hasattr(self.get_frame, 'clear_cache'):
            self.get_frame.clear_cache()

    def __iter__(self):
        return self.slice(self._slice)

    def __getitem__(self, item):
        if isinstance(item, slice):
            sliced = copy(self)
            sliced._slice = item
            return sliced
        else:
            return self.get_frame(item)

    def __len__(self):
        return len(self.frames)

    def __hash__(self):
        return merge_hashes(_hash(self.frames), _hash(self.atom_filter), _hash(self._slice))

    def __repr__(self):
        return "Coordinates <{}>: {}".format(self.frames.filename, self.atom_subset)

    @wraps(AtomSubset.subset)
    def subset(self, **kwargs):
        return Coordinates(self.frames, atom_subset=self.atom_subset.subset(**kwargs))

    @property
    def description(self):
        return self.atom_subset.description

    @description.setter
    def description(self, desc):
        self.atom_subset.description = desc


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
        self.frames = self.coordinates.frames
        self.atom_subset = self.coordinates.atom_subset
        self.function = function

    def __iter__(self):
        for frame in self.coordinates:
            yield self.function(frame)

    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.__class__(self.coordinates[item], self.function)
        else:
            frame = self.function(self.coordinates.__getitem__(item))
            if not isinstance(frame, CoordinateFrame):
                frame = frame.view(CoordinateFrame)
            frame.coordinates = self
            frame.step = item
            return frame

    def __len__(self):
        return len(self.coordinates.frames)

    def __hash__(self):
        return merge_hashes(_hash(self.coordinates), _hash(self.function))

    @wraps(Coordinates.subset)
    def subset(self, **kwargs):
        return CoordinatesMap(self.coordinates.subset(**kwargs), self.function)

    @property
    def description(self):
        return self.atom_subset.description

    @description.setter
    def description(self, desc):
        self.atom_subset.description = desc


class CoordinatesFilter:

    @property
    def atom_subset(self):
        pass

    def __init__(self, coordinates, atom_filter):
        self.coordinates = coordinates
        self.atom_filter = atom_filter

    def __getitem__(self, item):
        if isinstance(item, slice):
            sliced = copy(self)
            sliced.coordinates = self.coordinates[item]
            return sliced
        else:
            frame = self.coordinates[item]
            return frame[self.atom_filter]


class CoordinatesKDTree:
    """
    A KDTree of coordinates frames. The KDtrees are cached by a :func:`functools.lru_cache`.
    Uses :class:`scipy.spatial.cKDTree` by default, since it's significantly faster.
    Make sure to use scipy 0.17 or later or switch to the normal KDTree, since cKDTree has
    a memory leak in earlier versions.
    """

    def clear_cache(self):
        """Clear the LRU cache."""
        self._get_tree_at_index.cache_clear()

    @property
    def cache_info(self):
        """Return info about the state of the cache."""
        return self._get_tree_at_index.cache_info()

    def _get_tree_at_index(self, index):
        frame = self.frames[index]
        return self.kdtree(frame[self.selector(frame)])

    def __init__(self, frames, selector=None, maxcache=128, ckdtree=True):
        """
        Args:
            frames: Trajectory of the simulation, can be Coordinates object or reader
            selector: Selector function that selects a subset of each frame
            maxcache: Maxsize of the :func:`~functools.lru_cache`
            ckdtree: Use :class:`~scipy.spatial.cKDTree` or :class:`~scipy.spatial.KDTree` if False
        """
        if selector is not None:
            self.selector = selector
        else:
            self.selector = lambda x: slice(None)
        self.frames = frames
        self.kdtree = cKDTree if ckdtree else KDTree
        self._get_tree_at_index = lru_cache(maxsize=maxcache)(self._get_tree_at_index)

    def __getitem__(self, index):
        return self._get_tree_at_index(index)

    def __hash__(self):
        return merge_hashes(_hash(self.selector), _hash(self.frames))

    def __eq__(self, other):
        return super().__eq__(other)


def map_coordinates(func):
    @wraps(func)
    def wrapped(coordinates, *args, **kwargs):
        return CoordinatesMap(coordinates, partial(func, *args, **kwargs))
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
        origin: Origin of the pore which will be the coordinates origin after mapping
        sym_axis (opt.): Symmtery axis of the pore, may be a literal direction
            'x', 'y' or 'z' or an array of shape (3,)
    """
    if sym_axis in ('x', 'y', 'z'):
        rot_axis = np.zeros(shape=(3,))
        rot_axis[['x', 'y', 'z'].index(sym_axis)] = 1
    else:
        rot_axis = sym_axis

    return rotate_axis(coordinates - origin, rot_axis)


@map_coordinates
def vectors(coordinates, atoms_a, atoms_b, normed=False):
    """
    Compute the vectors between the atoms of two subsets.

    Args:
        coordinates: The Coordinates object the atoms will be taken from
        atoms_a: Mask or indices of the first atom subset
        atoms_b: Mask or indices of the second atom subset
        normed (opt.): If the vectors should be normed

    The defintion of atoms_a/b can be any possible subript of a numpy array.
    They can, for example, be given as a masking array of bool values with the
    same length as the frames of the coordinates. Or they can be a list of
    indices selecting the atoms of these indices from each frame.

    It is possible to compute the mean of several atoms before calculating the vectors,
    by using a two-dimensional list of indices. The following code computes the vectors
    between atoms 0, 3, 6 and the mean coordinate of atoms 1, 4, 7 and 2, 5, 8::

        >>> inds_a = [0, 3, 6]
        >>> inds_b = [[1, 4, 7], [2, 5, 8]]
        >>> vectors(coords, inds_a, inds_b)
        array([
            coords[0] - (coords[1] + coords[2])/2,
            coords[3] - (coords[4] + coords[5])/2,
            coords[6] - (coords[7] + coords[8])/2,
        ])
    """
    coords_a = coordinates[atoms_a]
    if len(coords_a.shape) > 2:
        coords_a = coords_a.mean(axis=0)
    coords_b = coordinates[atoms_b]
    if len(coords_b.shape) > 2:
        coords_b = coords_b.mean(axis=0)
    vectors = coords_a - coords_b
    norm = np.linalg.norm(vectors, axis=-1).reshape(-1, 1) if normed else 1
    return vectors / norm
