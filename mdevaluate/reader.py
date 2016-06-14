"""
Module that provides different readers for trajectory files.
"""

from .utils import hash_anything, merge_hashes

from functools import lru_cache

import pygmx
from pygmx.errors import InvalidMagicException, InvalidIndexException

def open(filename, cached=False):
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

    """

    if cached is not False:
        if isinstance(cached, bool):
            maxsize = 128
        else:
            maxsize = cached
        return CachedReader(filename, maxsize)
    else:
        return BaseReader(filename)


class BaseReader:
    """Base class for trajectory readers."""

    @property
    def filename(self):
        return self.rd.filename

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

    def __getitem__(self, item):
        return self.rd[item]

    def __len__(self):
        return len(self.rd)

    def __hash__(self):
        return merge_hashes(hash_anything(self.rd.filename),
                            hash_anything(str(self.rd.cache)))


class CachedReader(BaseReader):
    """A reader that has a least-recently-used cache for frames."""

    @property
    def cache_info(self):
        """Get Information about the lru cache."""
        return self._get_item.cache_info()

    def clear_cache(self):
        """Clear the cache of the frames."""
        self._get_item.cache_clear()

    def __init__(self, filename, maxsize):
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
