
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


def pbc_diff(v1, v2, box):
    """
    Calculate the difference of two vectors, considering periodic boundary conditions.
    """
    v = v1 - v2
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


def nojump(frame):
    """
    Return the nojump coordinates of a frame, based on a jump matrix.
    """
    selection = frame.coordinates.atom_subset.selection
    delta = np.array(np.vstack(
        [m[:frame.step + 1, selection].sum(axis=0) for m in frame.coordinates.frames.nojump_matrixes]
    ).T) * frame.box.diagonal()
    return frame - delta
