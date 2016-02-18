import numpy as np
from numpy import genfromtxt
import itertools
import re

from .reader import XTCReader, TRRReader
from .xtcindex import index_xtcfile

def atoms_from_grofile(file):
    """ Deprecated. Use Atoms.from_grofile"""
    t = genfromtxt(file, skip_header=2, delimiter=(10, 5), dtype='U8', autostrip=True, skip_footer=1)
    return t


def _atoms_from_grofile(file):
    t = genfromtxt(file, skip_header=2, delimiter=(5, 5, 5), dtype='U8', autostrip=True, skip_footer=1)
    return t


group_re = re.compile('\[ ([-+\w]+) \]')


def load_indices(indexfile):
    indices = {}
    index_array = None
    with open(indexfile) as idx_file:
        for line in idx_file:
            m = group_re.search(line)
            if m is not None:
                group_name = m.group(1)
                index_array = indices.get(group_name, [])
                indices[group_name] = index_array
            else:
                elements = line.strip().split('\t')
                elements = [x.split(' ') for x in elements]
                elements = itertools.chain(*elements)  # make a flat iterator
                elements = [x for x in elements if x != '']
                index_array += [int(x) - 1 for x in elements]
    return indices


def load_index(indexfile, index_names):
    """
    Load an indexfile, returns tuple of indices for index_names
    """
    index_arrays = [[] for _ in index_names]
    index_array = []  # This array will contain all indices before the first index_group. Should be empty

    with open(indexfile) as idx_file:
        for line in idx_file:
            if line.startswith("["):
                for i, index in enumerate(index_names):
                    if line.startswith("[ " + index + " ]"):
                        index_array = index_arrays[i]
                        break
                    index_array = None
            else:
                elements = line.strip().split('\t')
                elements = [x.split(' ') for x in elements]

                elements = itertools.chain(*elements)  # make a flat iterator
                elements = [x for x in elements if x != '']
                index_array += [int(x) - 1 for x in elements]
    return index_arrays
