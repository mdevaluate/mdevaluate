import numpy as np
from .atoms import AtomSubset
from .pbc import pbc_diff

class Coordinates:
    atom_filter = None

    def __init__(self, frames, atom_filter=None, atom_subset: AtomSubset=None):
        self.frames = frames
        assert atom_filter is None or atom_subset is None, "Cannot use both: subset and filter"

        if atom_filter is not None:
            self.atom_filter = atom_filter

        if atom_subset is not None:
            self.atom_filter = atom_subset.selection
            self.atom_subset = atom_subset
            self.atoms = atom_subset.atoms

    def __getitem__(self, item):
        try:
            if self.atom_filter is not None:
                return self.frames.__getitem__(item).coordinates[self.atom_filter]
            else:
                return self.frames.__getitem__(item).coordinates
        except EOFError:
            raise IndexError

    def __len__(self):
        return len(self.frames)


class MeanCoordinates(Coordinates):
    def __init__(self, frames, atom_filter=None, mean=1):
        super().__init__(frames, atom_filter)
        self.mean = mean
        assert mean >= 1, "Mean must be positive"

    def __getitem__(self, item):
        frame = super().__getitem__(item)
        for i in range(item+1, item+self.mean):
            frame += super().__getitem__(i)

        return frame/self.mean
    def len(self):
        return len(super()-self.mean+1)


class CoordinatesMap:
    def __init__(self, coordinates, function):
        self.coordinates = coordinates
        self.function = function
    
    def __getitem__(self, item):
        return self.function(self.coordinates.__getitem__(item))
    
    def __len__(self):
        return len(self.coordinates.frames)

def map_coordinates(fn):
    def wrapped(coordinates, *args, **kwargs):
        return CoordinatesMap(coordinates, lambda x:fn(x, *args, **kwargs))
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
    number_of_coordinates,number_of_dimensions = c.shape
    number_of_new_coordinates = number_of_coordinates // number_of_masses
    grouped_masses = c.reshape(number_of_new_coordinates, number_of_masses, number_of_dimensions)
    
    return np.average(grouped_masses, axis=1, weights=masses)
