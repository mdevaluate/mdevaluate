import os
from glob import glob

from . import atoms
from . import coordinates
from . import correlation
from . import distribution
from . import functions
from . import pbc
from . import simulation
from . import autosave
from . import reader

from .meta import notify

from pygmx import gromacs
from pygmx.errors import FileTypeError

__version__ = '16.09'


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


def open(directory, topology='*.tpr', trajectory='*.xtc',
         index_file=None, cached=False, reindex=True, verbose=True):
    """
    Open a simulation from a directory.

    Args:
        directory: Directory of the simulation.
        topology (opt.):
            Descriptor of the topology file (tpr or gro). By default a tpr file is
            used, if there is exactly one in the directoy.
        trajectory (opt.): Descriptor of the trajectory (xtc file).
        index_file (opt.): An index file that will be loaded alongside with the topology.
        cached (opt.):
            If the trajectory reader should be cached. Can be True, an integer or None.
            If this is True maxsize is 128, otherwise this is used as maxsize for
            the cache, None means infinite cache (this is a potential memory leak!).
        reindex (opt.): Regenerate the xtc-index if necessary.
        verbose (opt.): Be verbose about the opened files.

    Returns:
        A Coordinate object of the simulation.

    Example:
        Open a simulation located in '/path/to/sim', where the trajectory is
        located in a sub-directory '/path/to/sim/out' and named for Example
        'nojump_traj.xtc'. All read frames will be cached in memory.

        >>> open('/path/to/sim', trajectory='out/nojump*.xtc', cached=None)

    The file descriptors can use unix style pathname expansion to define the filenames.
    The default patterns use the recursive placeholder `**` which matches the base or
    any subdirctory, thus files in subdirectories with matching file type will be found too.
    For example: 'out/nojump*.xtc' would match xtc files in a subdirectory `out` that
    start with `nojump` and end with `.xtc`.

    For more details see: https://docs.python.org/3/library/glob.html
    """

    top_glob = glob(os.path.join(directory, topology), recursive=True)
    if top_glob is not None and len(top_glob) is 1:
        top_file, = top_glob
        top_ext = top_file.split('.')[-1]
        notify('Loading topology: {}'.format(top_file), verbose)
        if index_file is not None:
            index_glob = glob(os.path.join(directory, index_file), recursive=True)
            if index_glob is not None:
                index_file = index_glob[0]
            else:
                index_file = None
        if top_ext == 'tpr':
            atom_set = atoms.from_tprfile(top_file, index_file)
        elif top_ext == 'gro':
            atom_set = atoms.from_grofile(top_file, index_file)
        else:
            raise FileTypeError('Can not open file: {}'.format(top_file))

    else:
        raise FileNotFoundError('Topology file could not be identified.')

    traj_glob = glob(os.path.join(directory, trajectory), recursive=True)
    if traj_glob is not None and len(traj_glob) is 1:
        notify('Loading trajectory: {}'.format(traj_glob[0]), verbose)
        frames = reader.open(traj_glob[0], cached=cached, reindex=reindex)
    else:
        raise FileNotFoundError('Trajectory file could not be identified.')

    return coordinates.Coordinates(frames, atom_subset=atom_set)
