import os
from glob import glob

from . import atoms, coordinates, correlation, coordinates, distribution, functions, gromacs, pbc, simulation, autosave


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


def load_simulation(directory, xtc='*.xtc', tpr='*.tpr', gro='*.gro'):
    """
    Load a simulation from adirectory.
    """
    tpr_glob = glob(os.path.join(directory, tpr))
    gro_glob = glob(os.path.join(directory, gro))
    if tpr_glob is not None and len(tpr_glob) is 1:
        print('Loading topology: {}'.format(tpr_glob[0]))
        atom_set = atoms.from_tprfile(tpr_glob[0])
    elif gro_glob is not None and len(gro_glob) is 1:
        print('Loading topology: {}'.format(gro_glob[0]))
        atom_set = atoms.from_grofile(gro_glob[0])
    else:
        raise FileNotFoundError('Topology file could not be identified.')
    xtc_file, = glob(os.path.join(directory, xtc))
    print('Loading trajectory: {}'.format(xtc_file))
    frames = trajectory_from_xtc(xtc_file)

    return coordinates.Coordinates(frames, atom_subset=atom_set)
