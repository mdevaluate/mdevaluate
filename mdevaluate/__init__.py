import os
from glob import glob

from . import atoms
from . import coordinates
from . import correlation
from . import distribution
from . import functions
from . import pbc
from . import autosave
from . import reader
from .logging import logger

__version__ = '17.06'


def open(directory, topology='*.tpr', trajectory='*.xtc',
         index_file=None, cached=False, reindex=True, verbose=True, nojump=False,
         ignore_index_timestamps=False):
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
        nojump (opt.): If nojump matrixes should be generated. They will alwyas be loaded if present
        ignore_index_timestamps (opt.): 
            Ignore timestamps in xtc index file. If True the index file will be loaded 
            regardless of the timestamp

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
        logger.info('Loading topology: {}'.format(top_file))
        if index_file is not None:
            index_glob = glob(os.path.join(directory, index_file), recursive=True)
            if index_glob is not None:
                index_file = index_glob[0]
            else:
                index_file = None
    else:
        raise FileNotFoundError('Topology file could not be identified.')

    traj_glob = glob(os.path.join(directory, trajectory), recursive=True)
    if traj_glob is not None and len(traj_glob) is 1:
        traj_file = traj_glob[0]
        logger.info('Loading trajectory: {}'.format(traj_file))
    else:
        raise FileNotFoundError('Trajectory file could not be identified.')
    atom_set, frames = reader.open(
        top_file, traj_file, index_file=index_file, cached=cached, reindex=reindex,
        ignore_index_timestamps=ignore_index_timestamps
    )
    coords = coordinates.Coordinates(frames, atom_subset=atom_set)
    if nojump:
        try:
            frames.nojump_matrixes
        except reader.NojumpError:
            reader.generate_nojump_matrixes(coords)
    return coords


def open_energy(energyfile):
    """Open an energy file with EnergyReader."""
    edrfile = glob(energyfile)[0]
    return reader.EnergyReader(edrfile)
