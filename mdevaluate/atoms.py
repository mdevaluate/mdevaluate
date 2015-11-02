from scipy.spatial.distance import cdist
from .pbc import pbc_diff
from .gromacs import _atoms_from_grofile, load_indices
import numpy as np

from scipy.spatial import KDTree

def from_grofile(grofile, index_file=None):
    indices = None
    if index_file is not None:
        indices = load_indices(index_file)

    return Atoms(_atoms_from_grofile(grofile), indices).subset()

class Atoms:
    def __init__(self, atoms, indices=None):
        self.residue_ids, self.residue_names, self.atom_names = atoms.T
        self.residue_ids = np.array([int(m) for m in self.residue_ids])
        self.indices = indices

    def subset(self, *args, **kwargs):
        return AtomSubset(self).subset(*args, **kwargs)
    
    def __len__(self):
        return len(self.atom_names)

class AtomMismatch(Exception):
    pass

class AtomSubset:
    def __init__(self, atoms, selection=None):
        if selection is None:
            selection = np.ones(len(atoms), dtype='bool')
        self.selection = selection
        self.atoms = atoms

    def subset(self, atom_name=None, residue_name=None, residue_id=None, indices=None):
        new_subset = self
        if atom_name is not None:
            new_subset &= AtomSubset(self.atoms, self.atoms.atom_names == atom_name)

        if residue_name is not None:
            new_subset &= AtomSubset(self.atoms, self.atoms.residue_names == residue_name)

        if residue_id is not None:
            new_subset &= AtomSubset(self.atoms, self.atoms.residue_ids == residue_id)

        if indices is not None:
            selection = np.zeros(len(self.selection), dtype='bool')
            selection[indices] = True
            new_subset &= AtomSubset(self.atoms, selection)
        return new_subset

    @property
    def atom_names(self):
        return self.atoms.atom_names[self.selection]

    @property
    def residue_names(self):
        return self.atoms.residue_names[self.selection]

    @property
    def residue_ids(self):
        return self.atoms.residue_ids[self.selection]

    @property
    def indices(self):
        return np.where(self.selection)

    def __getitem__(self, slice):
        if isinstance(slice, str):
            indices = self.atoms.indices[slice]
            return self.atoms.subset()[indices] & self
        
        return self.subset(indices=self.indices[0].__getitem__(slice))

    def __and__(self, other):
        selection = (self.selection & other.selection)
        if self.atoms != other.atoms: raise AtomMismatch

        return AtomSubset(self.atoms, selection)

    def __or__(self, other):
        selection = (self.selection | other.selection)
        if self.atoms != other.atoms: raise AtomMismatch

        return AtomSubset(self.atoms, selection)

    def __invert__(self):
        selection = ~self.selection
        return AtomSubset(self.atoms, selection)

    def __repr__(self):
        return 'Subset of Atoms ({} of {})'.format(len(self.atoms.residue_names[self.selection]), len(self.atoms))
    
    @property
    def description(self):
        return "\n".join(["{}{} {}".format(resid, resname, atom_names) for resid, resname, atom_names in 
            zip(self.residue_ids, self.residue_names, self.atom_names)
        ])


def center_of_mass(position, mass=None):
    if mass is not None:
        return 1/mass.sum() * (mass*position).sum(axis=0)
    else:
        return 1/len(position) * position.sum(axis=0)


def gyration_radius(position, mass=None):
    r_s = center_of_mass(position, mass)
    
    return 1/len(position) * cdist(position, [r_s]).sum()
    

def layer_of_atoms(atoms, 
                   thickness,
                   plane_offset=np.array([0,0,0]),
                   plane_normal=np.array([1,0,0])):
                   
    p_ = atoms-plane_offset
    distance = np.dot(p_, plane_normal)
    
    return abs(distance) <= thickness

def distance_to_atoms(ref, atoms, box=None):
    """Get the minimal distance from atoms to ref.
    The result is an array of with length == len(atoms)
    """
    out = np.empty(atoms.shape[0])
    for i,atom in enumerate(atoms):
        diff = (pbc_diff(atom, ref, box) ** 2).sum(axis=1).min()
        out[i] = np.sqrt(diff)
    return out

def next_neighbors(atoms, number_of_neighbors=1, distance_upper_bound=np.inf):
    tree = KDTree(atoms)
    _, indices = tree.query(atoms, number_of_neighbors + 1, distance_upper_bound=distance_upper_bound)
    return indices[:,1:] # don't return the atoms itself