import os
from glob import glob

from . import atoms, coordinates, correlation, coordinates, distribution, evaluation, functions, gromacs, pbc, simulation, autosave


def trajectory_from_xtc(xtc_file, generate_index=True):
    """
    Load a trajectory from a xtc file. If no index file is found,
    a new one will be generated if generate_index is true.
    """
    index_file = gromacs.reader.index_filename_for_xtc(xtc_file)
    if not os.path.exists(index_file):
        print('No index file found, generating a new one. This may take a while...')
        gromacs.index_xtcfile(xtc_file)

    return gromacs.XTCReader(xtc_file)


def load_simulation(directory, xtc='*.xtc', gro='*.gro'):
    """
    Load a simulation from adirectory.
    """
    xtc_file, = glob(os.path.join(directory, xtc))
    gro_file, = glob(os.path.join(directory, gro))
    frames = trajectory_from_xtc(xtc_file)
    atom_set = atoms.from_grofile(gro_file)
    return coordinates.Coordinates(frames, atom_subset=atom_set)

#__all__ = ['atoms', 'coordinates']
