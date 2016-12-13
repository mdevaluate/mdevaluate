"""
Module that provides different readers for trajectory files.
"""

from .utils import hash_anything, merge_hashes

from functools import lru_cache
from collections import namedtuple
import os
from os import path
import logging
from array import array

import numpy as np
from scipy import sparse

import pygmx
from pygmx.errors import InvalidMagicException, InvalidIndexException


class NojumpError(Exception):
    pass


def open(filename, cached=False, reindex=False):
    """
    Opens a trajectory file with the apropiate reader.

    Args:
        filename (str):
            Trajectory file to open, the reader will be chosen
            according to the file extension.
        cached (opt.):
            If Reader should be cached with lru_cache. If this is True, maxsize for
            the cache is 128, otherwise the argument is passed as maxsize.
            Use cached=None to get an unbound cache.
        reindex (opt.): Regenerate the index of the xtc-file
        nojump (opt.): If nojump matrixes should be generated.

    """

    if cached is not False:
        if isinstance(cached, bool):
            maxsize = 128
        else:
            maxsize = cached
        reader = CachedReader(filename, maxsize, reindex=reindex)
    else:
        reader = BaseReader(filename, reindex=reindex)

    return reader


def nojump_filename(reader):
    directory, fname = path.split(reader.filename)
    return path.join(directory, '.{}.nojump.npz'.format(fname))

CSR_ATTRS = ('data', 'indices', 'indptr')
NOJUMP_MAGIC = 2016


def parse_jumps(trajectory):
    prev = trajectory[0].whole
    box = prev.box.diagonal()
    SparseData = namedtuple('SparseData', ['data', 'row', 'col'])
    jump_data = (
        SparseData(data=array('b'), row=array('l'), col=array('l')),
        SparseData(data=array('b'), row=array('l'), col=array('l')),
        SparseData(data=array('b'), row=array('l'), col=array('l'))
    )

    for i, curr in enumerate(trajectory):
        if i % 500 == 0:
            logging.debug('Parse jumps Step: %d', i)
        delta = ((curr - prev) / box).round().astype(np.int8)
        prev = curr
        for d in range(3):
            col, = np.where(delta[:, d] != 0)
            jump_data[d].col.extend(col)
            jump_data[d].row.extend([i] * len(col))
            jump_data[d].data.extend(delta[col, d])

    return jump_data


def generate_nojump_matrixes(trajectory):
    """
    Create the matrixes with pbc jumps for a trajectory.
    """
    logging.info('generate Nojump Matrixes for: {}'.format(trajectory))

    jump_data = parse_jumps(trajectory)
    N = len(trajectory)
    M = len(trajectory[0])

    trajectory.frames.nojump_matrixes = tuple(
        sparse.csr_matrix((np.array(m.data), (m.row, m.col)), shape=(N, M)) for m in jump_data
    )


def save_nojump_matrixes(reader, matrixes=None):
    if matrixes is None:
        matrixes = reader.nojump_matrixes
    data = {'checksum': merge_hashes(NOJUMP_MAGIC, hash(reader))}
    for d, mat in enumerate(matrixes):
        data['shape'] = mat.shape
        for attr in CSR_ATTRS:
            data['{}_{}'.format(attr, d)] = getattr(mat, attr)

    np.savez(nojump_filename(reader), **data)


def load_nojump_matrixes(reader):
    try:
        data = np.load(nojump_filename(reader))
    except AttributeError:
        # npz-files can be corrupted, propably a bug for big arrays saved with savez_compressed?
        os.remove(nojump_filename(reader))
        return
    if data['checksum'] == merge_hashes(NOJUMP_MAGIC, hash(reader)):
        reader.nojump_matrixes = tuple(
            sparse.csr_matrix(
                tuple(data['{}_{}'.format(attr, d)] for attr in CSR_ATTRS),
                shape=data['shape']
            )
            for d in range(3)
        )
        logging.info('Loaded Nojump Matrixes: {}'.format(nojump_filename(reader)))
    else:
        logging.info('Invlaid Nojump Data: {}'.format(nojump_filename(reader)))


class BaseReader:
    """Base class for trajectory readers."""

    @property
    def filename(self):
        return self.rd.filename

    @property
    def nojump_matrixes(self):
        if self._nojump_matrixes is None:
            raise NojumpError('Nojump Data not available: {}'.format(self.filename))
        return self._nojump_matrixes

    @nojump_matrixes.setter
    def nojump_matrixes(self, mats):
        self._nojump_matrixes = mats
        save_nojump_matrixes(self)

    def __init__(self, filename, reindex=False):
        """
        Args:
            filename: Trajectory file to open.
            reindex (bool, opt.): If True, regenerate the index file if necessary.
        """
        try:
            self.rd = pygmx.open(filename)
        except InvalidMagicException:
            raise InvalidIndexException('This is not a valid index file: {}'.format(filename))
        except InvalidIndexException:
            if reindex:
                pygmx.gromacs.index_xtcfile(filename)
                self.rd = pygmx.open(filename)
            else:
                raise InvalidIndexException('Index file is invalid, us reindex=True to regenerate.')

        self._nojump_matrixes = None
        if path.exists(nojump_filename(self)):
            load_nojump_matrixes(self)

    def __getitem__(self, item):
        return self.rd[item]

    def __len__(self):
        return len(self.rd)

    def __hash__(self):
        hashes = [hash_anything(self.rd.filename), hash_anything(str(self.rd.cache))]
        return merge_hashes(*hashes)


class CachedReader(BaseReader):
    """A reader that has a least-recently-used cache for frames."""

    @property
    def cache_info(self):
        """Get Information about the lru cache."""
        return self._get_item.cache_info()

    def clear_cache(self):
        """Clear the cache of the frames."""
        self._get_item.cache_clear()

    def __init__(self, filename, maxsize, reindex=False):
        """
        Args:
            filename (str): Trajectory file that will be opened.
            maxsize: Maximum size of the lru_cache or None for infinite cache.
        """
        super().__init__(filename)
        self._get_item = lru_cache(maxsize=maxsize)(self._get_item)

    def _get_item(self, item):
        """Buffer function for lru_cache, since __getitem__ can not be cached."""
        return super().__getitem__(item)

    def __getitem__(self, item):
        return self._get_item(item)


class EnergyReader:
    """A reader for Gromacs energy files."""

    def __init__(self, edrfile):
        """
        Args:
            edrfile: Filename of the energy file
            topology (opt.): Filename of the topology, speeds up file io since the length of the energy file is known
        """
        edr = pygmx.open(edrfile)
        self.time, data = edr.read()
        self.types, self.units = zip(*edr.types)
        self.data = data.T

    def __getitem__(self, type):
        """
        Get time series of an energy type.
        """
        if type in self.types:
            return self.data[self.types.index(type)]
        else:
            raise KeyError('Energy type {} not found in Energy File.'.format(type))
