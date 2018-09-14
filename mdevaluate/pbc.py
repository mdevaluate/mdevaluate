
from collections import OrderedDict

import numpy as np
import numba

from scipy.spatial import cKDTree

from .logging import logger

def pbc_diff_old(v1, v2, box):
    """
    Calculate the difference of two vestors, considering optional boundary conditions.
    """
    if box is None:
        v = v1 - v2
    else:
        v = v1 % box - v2 % box
        v -= (v > box / 2) * box
        v += (v < -box / 2) * box

    return v


def pbc_diff(v1, v2=None, box=None):
    """
    Calculate the difference of two vectors, considering periodic boundary conditions.
    """
    if v2 is None:
        v = v1
    else:
        v = v1 -v2
    if box is not None:
        s = v / box
        v = box * (s - s.round())
    return v


@numba.jit(nopython=True)
def pbc_diff_numba(ri, rj, box):
    v = ri % box - rj % box
    v -= (v > box / 2) * box
    v += (v < -box / 2) * box
    return v


def whole(frame, residue_ids=None, len_res=None):
    """
    Apply ``-pbc whole`` to a CoordinateFrame.
    At the moment this function needs additonal information about the bonds within the residue.

    Args:
        frame: The original CoordinateFrame
        bonds:
            Information about bonds within the residue, as a list of shape (N, 3) where
            each entry of the list defines tweo bonded atoms and the bond length.
        residue_ids: The residue_ids of the atoms, if None this is taken from the CoordinateFrame

    """
    if residue_ids is None:
        residue_ids = frame.coordinates.atom_subset.residue_ids
    if len_res is None:
        len_res = (residue_ids == residue_ids[0]).sum()
    box = frame.box.diagonal()

    nr_res = len(frame) // len_res
    residues = frame.reshape(nr_res, len_res, 3)
    masses = frame.masses.reshape(nr_res, len_res, 1) / frame.masses[:len_res].sum()
    com = (residues * masses).sum(axis=1).reshape(-1, 1, 3)
    com_dist = residues - com

    correction = np.zeros(residues.shape)
    n, m, d = np.where(com_dist > box / 2)
    correction[n, m, d] = -box[d]
    n, m, d = np.where(com_dist < -box / 2)
    correction[n, m, d] = box[d]

    if len_res == 2 and frame.masses[0] == frame.masses[1]:
        # special case for 2-site molecules, where the algorithm fails since the COM is somewhere in the box center
        correction[:,1,:] = 0

    whole_frame = residues + correction
    return whole_frame.reshape(nr_res * len_res, 3)


NOJUMP_CACHESIZE = 128


def nojump(frame, usecache=True):
    """
    Return the nojump coordinates of a frame, based on a jump matrix.
    """
    selection = frame.coordinates.atom_subset.selection

    reader = frame.coordinates.frames
    if usecache:
        if not hasattr(reader, '_nojump_cache'):
            reader._nojump_cache = OrderedDict()
        i0s = [x for x in reader._nojump_cache if x <= frame.step]
        if len(i0s) > 0:
            i0 = max(i0s)
            delta = reader._nojump_cache[i0]
            i0 += 1
        else:
            i0 = 0
            delta = 0

        delta += np.array(np.vstack(
            [m[i0:frame.step + 1].sum(axis=0) for m in frame.coordinates.frames.nojump_matrixes]
        ).T) * frame.box.diagonal()

        reader._nojump_cache[frame.step] = delta
        while len(reader._nojump_cache) > NOJUMP_CACHESIZE:
            reader._nojump_cache.popitem(last=False)
        delta = delta[selection, :]
    else:
        delta = np.array(np.vstack(
            [m[:frame.step + 1, selection].sum(axis=0) for m in frame.coordinates.frames.nojump_matrixes]
        ).T) * frame.box.diagonal()
    return frame - delta


def pbc_points(coordinates, box, thickness=0, index=False, inclusive=True, center=None):
    """
    Returns the points their first periodic images. Does not fold them back into the box.
    Thickness 0 means all 27 boxes. Positive means the box+thickness. Negative values mean that less than the box is returned.
    index=True also returns the indices with indices of images being their originals values.
    inclusive=False returns only images, does not work with thickness <= 0
    """
    if center is None:
        center = box/2
    allcoordinates = np.copy(coordinates)
    indices = np.tile(np.arange(len(coordinates)),(27))
    for x in range(-1, 2, 1):
            for y in range(-1, 2, 1):
                for z in range(-1, 2, 1):
                    vv = np.array([x, y, z], dtype=float)
                    if not (vv == 0).all() :
                        allcoordinates = np.concatenate((allcoordinates, coordinates + vv*box), axis=0)
    
    if thickness != 0:
        mask = np.all(allcoordinates < center+box/2+thickness, axis=1)
        allcoordinates = allcoordinates[mask]
        indices = indices[mask]
        mask = np.all(allcoordinates > center-box/2-thickness, axis=1)
        allcoordinates = allcoordinates[mask]
        indices = indices[mask]
    if not inclusive and thickness > 0:
        allcoordinates = allcoordinates[len(coordinates):]
        indices = indices[len(coordinates):]
    if index:
        return (allcoordinates, indices)
    return allcoordinates
