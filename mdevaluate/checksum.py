
import functools
import hashlib
from .logging import logger

import numpy as np

# This variable is used within the checksum function to salt the md5 sum.
# May be changed to force a different checksum for similar objects.
SALT = 42


def version(version_nr, calls=[]):
    """Function decorator that assigns a custom checksum to a function."""

    def decorator(func):
        cs = checksum(func.__name__, version_nr, *calls)
        func.__checksum__ = lambda: cs

        @functools.wraps(func)
        def wrapped(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapped
    return decorator


def checksum(*args):
    """
    Calculate the checksum of any object.
    """
    bstr = str(SALT).encode()
    for arg in args:
        if hasattr(arg, '__checksum__'):
            logger.debug('Checksum via __checksum__: %s', str(arg))
            bstr += str(arg.__checksum__()).encode()
        elif arg is None:
            bstr += b'None'
        elif isinstance(arg, bytes):
            bstr += arg
        elif isinstance(arg, str):
            bstr += arg.encode()
        elif hasattr(arg, '__code__'):
            logger.debug('Checksum via __code__ for %s', str(arg))
            bstr += arg.__code__.co_code
            if arg.__closure__ is not None:
                for cell in arg.__closure__:
                    bstr += str(checksum(cell.cell_contents)).encode()
        elif isinstance(arg, functools.partial):
            logger.debug('Checksum via partial for %s', str(arg))
            bstr += str(checksum(arg.func)).encode()
            for x in arg.args:
                bstr += str(checksum(x)).encode()
            for k in sorted(arg.keywords.keys()):
                bstr += k.encode() + str(checksum(arg.keywords[k])).encode()
        elif isinstance(arg, np.ndarray):
            bstr += arg.tobytes()
        else:
            logger.debug('Checksum via str for %s', str(arg))
            bstr += str(arg).encode()

    m = hashlib.md5()
    m.update(bstr)
    return int.from_bytes(m.digest(), 'big')
