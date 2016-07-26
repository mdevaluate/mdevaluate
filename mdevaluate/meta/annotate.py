import warnings
import functools


def unfinished(f, comment=""):

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        warnings.warn('Using incomplete function {}: '.format(f.__name__, comment))
        return f(*args, **kwargs)
    return wrapped


def untested(f, comment=""):

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        warnings.warn('Using untested function {}: '.format(f.__name__, comment))
        return f(*args, **kwargs)
    return wrapped


def deprecated(replacement):
    def decorator(f):
        msg = '{} is deprecated, use {}.{} instead.'.format(f.__name__, replacement.__module, replacement.__name__)

        @functools.wraps(f)
        def wrapped(*args, **kwargs):
            warnings.warn(msg, DeprecationWarning, stacklevel=2)
            return f(*args, **kwargs)

        return wrapped
