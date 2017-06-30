
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


@numba.jit(nopython=True)
def pbc_diff_numba(ri, rj, box):
    v = ri % box - rj % box
    v -= (v > box / 2) * box
    v += (v < -box / 2) * box
    return v

def pbc_vec(a, b, box):
    """
    Calculate the vector pointing from a to b. This assumes that periodic
    boundary conditions apply and vectors are never longer than half the
    box size in any direction. Instead they will point into negative directions.
    """
    d = a-b
    d = d-(d//box)*box
    d = d-(d//(box/2))*box
    return d

def pbc_diff(a, b, box=None):
    """
    Calculate the difference of two vestors, considering optional boundary conditions.
    """
    if box is None:
        return a-b
    d = a-b
    d = d-(d//box)*box
    d = d-(d//(box/2))*box
    return d

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

    whole_frame = residues + correction
    return whole_frame.reshape(nr_res * len_res, 3)


def nojump(frame):
    """
    Return the nojump coordinates of a frame, based on a jump matrix.
    """
    selection = frame.coordinates.atom_subset.selection
    delta = np.array(np.vstack(
        [m[:frame.step + 1, selection].sum(axis=0) for m in frame.coordinates.frames.nojump_matrixes]
    ).T) * frame.box.diagonal()
    return frame - delta


def pbc_points(points, box, thickness=0, index=False, inclusive=True):
    """
    Returns the points folded back into the box and the first periodic images.
    Thickness 0 means all 27 boxes. Positive means the box+thickness. Negative values mean that less than the box is returned.
    index=True also returns the indices with indices of images being their originals values.
    inclusive=False returns only images, does not work with thickness <= 0
    """
    coordinates = np.copy(points)%box
    allcoordinates = np.copy(coordinates)
    indices = np.tile(np.arange(len(points)),(27))
    for x in range(-1, 2, 1):
            for y in range(-1, 2, 1):
                for z in range(-1, 2, 1):
                    vv = np.array([x, y, z], dtype=float)
                    if not (vv == 0).all() :
                        allcoordinates = np.concatenate((allcoordinates, coordinates + vv*box), axis=0)
    
    if thickness != 0:
        mask = np.all(allcoordinates < box+thickness, axis=1)
        allcoordinates = allcoordinates[mask]
        indices = indices[mask]
        mask = np.all(allcoordinates > -thickness, axis=1)
        allcoordinates = allcoordinates[mask]
        indices = indices[mask]
    if not inclusive and thickness > 0:
        allcoordinates = allcoordinates[len(points):]
        indices = indices[len(points):]
    if index:
        return (allcoordinates, indices)
    return allcoordinates
