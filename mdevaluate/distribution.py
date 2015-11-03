import numpy as np
from .coordinates import pbc_diff
from .atoms import next_neighbors


def time_average(function, coordinates, pool=None):

    coordinate_iter = iter(coordinates)
    number_of_averages = 0
    result = 0

    if pool is not None:
        _map = pool.imap
    else:
        _map = map

    evaluated = _map(function, coordinate_iter)

    for ev in evaluated:
        number_of_averages += 1
        result += ev
        if number_of_averages % 100 == 0:
            print(number_of_averages)

    return result/number_of_averages
    
def time_histogram(function, coordinates, bins, hist_range, pool=None):
    coordinate_iter = iter(coordinates)
    
    if pool is not None:
        _map = pool.imap
    else:
        _map = map

    evaluated = _map(function, coordinate_iter)
    
    results = []
    hist_results = []
    for num, ev in enumerate(evaluated):
        results.append(ev)
        
        if num % 100 == 0 and num > 0:
            print(num)
            r = np.array(results).T
            for i, row in enumerate(r):
                histo, _ = np.histogram(row, bins=bins, range=hist_range)
                if len(hist_results) <= i:
                    hist_results.append(histo)
                else:
                    hist_results[i] += histo
                results = []
    return hist_results


def radial_pair_distribution(atoms, bins, box=None):
    number_of_atoms = len(atoms)
    # Calculate the upper triangle of differences
    indices = np.triu_indices(number_of_atoms, k=1)

    # vector between atoms
    diff = pbc_diff(atoms[indices[0]], atoms[indices[1]], box)

    dist = (diff**2).sum(axis=1)**.5
    volume = 4/3*np.pi*(bins[1:]**3-bins[:-1]**3)
    hist = np.histogram(dist, bins)[0]
    density = (number_of_atoms-1) / np.prod(box)

    return hist / volume / density * (2/(number_of_atoms-1))


def distance_distribution(atoms, bins):
    connection_vectors = atoms[:-1,:]-atoms[1:,:]
    connection_lengths = (connection_vectors**2).sum(axis=1)**.5
    return np.histogram(connection_lengths, bins)[0]

def tetrahedral_order(atoms):
    indices = next_neighbors(atoms, number_of_neighbors = 4)
    neighbors = atoms[indices]
    neighbors_1, neighbors_2, neighbors_3, neighbors_4 = \
        neighbors[:,0,:], neighbors[:,1,:], neighbors[:,2,:], neighbors[:,3,:]
    
    # Connection vectors
    neighbors_1 -= atoms
    neighbors_2 -= atoms
    neighbors_3 -= atoms
    neighbors_4 -= atoms
    
    # Normed Connection vectors
    neighbors_1 /= ((neighbors_1**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0],1)
    neighbors_2 /= ((neighbors_2**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0],1)
    neighbors_3 /= ((neighbors_3**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0],1)
    neighbors_4 /= ((neighbors_4**2).sum(axis=1)**.5).reshape(neighbors_1.shape[0],1)
    
    a_1_2 = ((neighbors_1*neighbors_2).sum(axis=1) + 1/3)**2
    a_1_3 = ((neighbors_1*neighbors_3).sum(axis=1) + 1/3)**2
    a_1_4 = ((neighbors_1*neighbors_4).sum(axis=1) + 1/3)**2
    
    a_2_3 = ((neighbors_2*neighbors_3).sum(axis=1) + 1/3)**2
    a_2_4 = ((neighbors_2*neighbors_4).sum(axis=1) + 1/3)**2
     
    a_3_4 = ((neighbors_3*neighbors_4).sum(axis=1) + 1/3)**2
    
    q = 1 - 3/8 * (a_1_2 + a_1_3 + a_1_4 + a_2_3 + a_2_4 + a_3_4)
    
    return q
