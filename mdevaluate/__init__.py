import os
from glob import glob

from . import atoms, coordinates, correlation, distribution, functions, pbc, simulation, autosave

from pygmx import gromacs

__version__ = '1.2'

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


def load_simulation(directory, xtc='*.xtc', tpr='*.tpr', gro='*.gro', index_file=None, caching=False):
    """
    Load a simulation from a directory.

    Args:
        directory:  The directory where the simulation is located
        xtc (opt.): Descriptor of the trajectory file
        tpr (opt.): Descriptors of the tpr file, if this is None no tpr file will be used.
        gro (opt.): Descriptors of the gro file, if this is None no gro file will be used.
        index_file: Name of a index file that will be loaded with the gro file
        caching:    If caching should be activated in the Coordinates

    Only one topology filename has to be specified, tpr files will be prefered.
    The function uses :func:`trajectory_from_xtc` to load the xtc file, hence a new
    xtc-index file will be generated if necessary.

    The file descriptors can use unix style pathname expansion to define the filenames.
    For example: 'out/nojump*.xtc' would match xtc files in a subdirectory `out` that
    start with `nojump` and end with `.xtc`.

    For more details see: https://docs.python.org/3/library/glob.html
    """

    tpr_glob = glob(os.path.join(directory, tpr)) if tpr is not None else None
    gro_glob = glob(os.path.join(directory, gro)) if gro is not None else None
    if tpr_glob is not None and len(tpr_glob) is 1:
        print('Loading topology: {}'.format(tpr_glob[0]))
        atom_set = atoms.from_tprfile(tpr_glob[0])
    elif gro_glob is not None and len(gro_glob) is 1:
        print('Loading topology: {}'.format(gro_glob[0]))
        atom_set = atoms.from_grofile(gro_glob[0], index_file)
    else:
        raise FileNotFoundError('Topology file could not be identified.')
    xtc_file, = glob(os.path.join(directory, xtc))
    print('Loading trajectory: {}'.format(xtc_file))
    frames = trajectory_from_xtc(xtc_file)

    return coordinates.Coordinates(frames, atom_subset=atom_set, caching=caching)
