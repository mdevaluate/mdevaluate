"""
Collection of utility functions.
"""

def hash_of_iterable(iterable):
    """Return a hash value for the current state of an iterable."""
    return hash(''.join(str(x) for x in iterable))

def merge_hashes(*hashes):
    """Merge several hashes to one hash value."""
    return hash_of_iterable(hashes)
