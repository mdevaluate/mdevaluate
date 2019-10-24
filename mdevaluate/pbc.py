
from collections import OrderedDict
import os

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
    if box is None:
        out = v1 - v2
    elif len(getattr(box, 'shape', [])) == 1: 
        out = pbc_diff_rect(v1, v2, box)
    elif len(getattr(box, 'shape', [])) == 2: 
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


# Parameter used to switch reference position for whole molecules
# 'com': use center of mass of each molecule
# 'gmx': use Gromacs internal procedure
# 'simple': use first atom in each molecule 
WHOLEMODE = 'com'

fname = os.path.expanduser('~/.mdevaluate/WHOLEMODE')
if os.path.exists(fname):
    with open(fname) as f:
        WHOLEMODE = f.read().strip()
    logger.info('Setting WHOLEMODE to %s, according to file ~/.mdevaluate/WHOLEMODE', WHOLEMODE)

def whole(frame):
    """
    Apply ``-pbc whole`` to a CoordinateFrame.
    """
    residue_ids = frame.coordinates.atom_subset.residue_ids
    box = frame.box.diagonal()

    if WHOLEMODE == 'com':
        logger.debug('Using COM as reference for whole.')
        coms = np.array([
            np.bincount(residue_ids, weights=c * frame.masses)[1:] / np.bincount(residue_ids, weights=frame.masses)[1:]
            for c in frame.T
        ]).T[residue_ids - 1]
    elif WHOLEMODE == 'gmx':
        import pygmx
        
        x = pygmx.make_xtcframe_whole(
            frame.coordinates.frames[frame.step].positions,
            frame.box,
            frame.coordinates.atoms.reader
        )
        cor = x[frame.selection] - frame

    else:
        # make sure, residue_ids are sorted, then determine indices at which the res id changes
        # kind='stable' assures that any existent ordering is preserved
        logger.debug('Using first atom as reference for whole.')
        sort_ind = residue_ids.argsort(kind='stable')
        i = np.concatenate([[0], np.where(np.diff(residue_ids[sort_ind]) > 0)[0] + 1])
        coms = frame[sort_ind[i]][residue_ids - 1]
    
    if WHOLEMODE != 'gmx':
        cor = np.zeros_like(frame)
        cd = frame - coms
        n, d = np.where(cd > box / 2 * 0.9)
        cor[n, d] = -box[d]
        n, d = np.where(cd < -box / 2 * 0.9)
        cor[n, d] = box[d]

    # this fix is only necessary when COM is the reference
    if WHOLEMODE == 'com':
        duomask = np.bincount(residue_ids)[1:][residue_ids - 1] == 2
        if np.any(duomask):
            duomask[::2] = False
            cor[duomask] = 0

    return frame + cor


NOJUMP_CACHESIZE = 128


def nojump(frame, usecache=True):
    """
    Return the nojump coordinates of a frame, based on a jump matrix.
    """
    selection = frame.selection
    reader = frame.coordinates.frames
    if usecache:
        if not hasattr(reader, '_nojump_cache'):
            reader._nojump_cache = OrderedDict()
        # make sure to use absolute (non negative) index
        abstep = frame.step % len(frame.coordinates)
        i0s = [x for x in reader._nojump_cache if x <= abstep]
        if len(i0s) > 0:
            i0 = max(i0s)
            delta = reader._nojump_cache[i0]
            i0 += 1
        else:
            i0 = 0
            delta = 0

        delta = delta + np.array(np.vstack(
            [m[i0:abstep + 1].sum(axis=0) for m in reader.nojump_matrixes]
        ).T) * frame.box.diagonal()

        reader._nojump_cache[abstep] = delta
        while len(reader._nojump_cache) > NOJUMP_CACHESIZE:
            reader._nojump_cache.popitem(last=False)
        delta = delta[selection, :]
    else:
        delta = np.array(np.vstack(
            [m[:frame.step + 1, selection].sum(axis=0) for m in reader.nojump_matrixes]
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
