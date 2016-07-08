
import numpy as np


def pbc_diff(v1, v2, box):
    v = v1-v2
    if box is None:
        return v

    v -= (v > box/2)*box
    v += (v < -box/2) * box

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
    if len(np.unique(box)) == 1:
        box = box[0]
        correction[com_dist > box/2] = -box
        correction[com_dist < -box/2] = box
    else:
        raise NotImplementedError('Non cubic box is not implemented.')

    whole_frame = residues + correction
    return whole_frame.reshape(nr_res*len_res, 3)
