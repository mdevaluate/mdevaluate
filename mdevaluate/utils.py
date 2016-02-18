"""
Collection of utility functions.
"""


def hash_anything(arg):
    """Return a hash value for the current state of any argument."""
    try:
        return hash(arg)
    except TypeError:
        return hash(str(arg))


def merge_hashes(*hashes):
    """Merge several hashes to one hash value."""
    return hash_of_iterable(hashes)
