__author__ = 'mbartelm'
import logging
import functools


def unfinished(f, comment=""):

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        logging.warning('Using incomplete function {}: '.format(f.__name__, comment))
        return f(*args, **kwargs)
    return wrapped


def untested(f, comment=""):

    @functools.wraps(f)
    def wrapped(*args, **kwargs):
        logging.warning('Using untested function {}: '.format(f.__name__, comment))
        return f(*args, **kwargs)
    return wrapped
