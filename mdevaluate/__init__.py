from . import atoms, coordinates, correlation, coordinates, distribution, evaluation, functions, gromacs, pbc, simulation, autosave
import os

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


__all__ = ['atoms', 'coordinates']
