
import numpy as np
import numba

from scipy.spatial import cKDTree
from itertools import product

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
    if box==None:
        out = v1 - v2
    elif hasattr(box, 'shape') and len(box.shape) == 1: 
        out = pbc_diff_rect(v1, v2, box)
    elif hasattr(box, 'shape') and len(box.shape) == 2: 
        out = pbc_diff_tric(v1, v2, box)        
    else: raise NotImplementedError("cannot handle box")
    return out

def pbc_diff_rect(v1, v2, box):
    """
    Calculate the difference of two vectors, considering periodic boundary conditions.
    """
    if v2 is None:
        v = v1
    else:
        v = v1 -v2
   
    s = v / box
    v = box * (s - s.round())
    return v


def pbc_diff_tric(v1, v2=None, box=None):
    """
    difference vector for arbitrary pbc
    
        Args:
        box_matrix: CoordinateFrame.box
    """
    if len(box.shape) == 1: box = np.diag(box)  
    if v1.shape == (3,): v1 = v1.reshape((1,3)) #quick 'n dirty
    if v2.shape == (3,): v2 = v2.reshape((1,3))   
    if box is not None:
        r3 = np.subtract(v1, v2)
        r2 = np.subtract(r3, (np.rint(np.divide(r3[:,2],box[2][2])))[:,np.newaxis] * box[2][np.newaxis,:])
        r1 = np.subtract(r2, (np.rint(np.divide(r2[:,1],box[1][1])))[:,np.newaxis] * box[1][np.newaxis,:])
        v  = np.subtract(r1, (np.rint(np.divide(r1[:,0],box[0][0])))[:,np.newaxis] * box[0][np.newaxis,:])
    else:
        v = v1 - v2
    return v


def pbc_dist(a1,a2,box = None):
    return ((pbc_diff(a1,a2,box)**2).sum(axis=1))**0.5


def pbc_extend(c, box):
    """
     in: c is frame, box is frame.box
     out: all atoms in frame and their perio. image (shape => array(len(c)*27,3))
    """
    c=np.asarray(c)
    if c.shape == (3,): c = c.reshape((1,3)) #quick 'n dirty
    comb = np.array([np.asarray(i) for i in product([0,-1,1],[0,-1,1],[0,-1,1])])
    b_matrices = comb[:,:,np.newaxis]*box[np.newaxis,:,:]
    b_vectors = b_matrices.sum(axis=1)[np.newaxis,:,:]
    return (c[:,np.newaxis,:]+b_vectors)
    

def pbc_kdtree(v1,box, leafsize = 32, compact_nodes = False, balanced_tree = False):
    """
    kd_tree with periodic images
    box - whole matrix
    rest: optional optimization
    """
    r0 =  cKDTree(pbc_extend(v1,box).reshape((-1,3)),leafsize ,compact_nodes ,balanced_tree)
    return r0


def pbc_kdtree_query(v1,v2,box,n = 1):
    """
    kd_tree query with periodic images
    """
    r0, r1 =  pbc_kdtree(v1,box).query(v2,n)
    r1 = r1 // 27
    return r0, r1


def pbc_backfold_rect(act_frame,box_matrix):
    """
    mimics "trjconv ... -pbc atom -ur rect" 
    
    folds coords of act_frame in cuboid 
    
    """
    af=np.asarray(act_frame)
    if af.shape == (3,): act_frame = act_frame.reshape((1,3)) #quick 'n dirty
    b = box_matrix
    c = np.diag(b)/2
    af = pbc_diff(np.zeros((1,3)),af-c,b)
    return af + c


def pbc_backfold_compact(act_frame,box_matrix):
    """
    mimics "trjconv ... -pbc atom -ur compact" 
    
    folds coords of act_frame in wigner-seitz-cell (e.g. dodecahedron)
    """
    c = act_frame
    box = box_matrix
    ctr = box.sum(0)/2
    c=np.asarray(c)
    shape = c.shape
    if shape == (3,): 
        c = c.reshape((1,3))
        shape = (1,3)  #quick 'n dirty
    comb = np.array([np.asarray(i) for i in product([0,-1,1],[0,-1,1],[0,-1,1])])
    b_matrices = comb[:,:,np.newaxis]*box[np.newaxis,:,:]
    b_vectors = b_matrices.sum(axis=1)[np.newaxis,:,:]
    sc = c[:,np.newaxis,:]+b_vectors
    w = np.argsort((((sc)-ctr)**2).sum(2),1)[:,0]
    return sc[range(shape[0]),w]


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
