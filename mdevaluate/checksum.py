
import functools
import hashlib
from .logging import logger
from types import ModuleType, FunctionType
import inspect

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


def strip_comments(s):
    """Strips comment lines and docstring from Python source string."""
    o = ''
    in_docstring = False
    for l in s.split('\n'):
        if l.strip().startswith(('#', '"', "'")) or in_docstring:
            in_docstring = l.strip().startswith(('"""', "'''")) + in_docstring == 1
            continue
        o += l + '\n'
    return o


def checksum(*args, csum=None):
    """
    Calculate a checksum of any object, by md5 hash.

    Input for the hash are some salt bytes and the byte encoding of a string 
    that depends on the object and its type:

    - If a method __checksum__ is available, it's return value is converted to bytes
    - str or bytes are used as md5 input directly
    - modules use the __name__ attribute
    - functions use the function code and any closures of the function
    - functools.partial uses the checksum of the function and any arguments, that were defined
    - numpy.ndarray uses bytes representation of the array (arr.tobytes())
    - Anything else is converted to a str
    """
    if csum is None:
        csum = hashlib.sha1()
        csum.update(str(SALT).encode())

    for arg in args:
        if hasattr(arg, '__checksum__'):
            logger.debug('Checksum via __checksum__: %s', str(arg))
            csum.update(str(arg.__checksum__()).encode())
        elif isinstance(arg, bytes):
            csum.update(arg)
        elif isinstance(arg, str):
            csum.update(arg.encode())
        elif isinstance(arg, ModuleType):
            csum.update(arg.__name__.encode())
        elif isinstance(arg, FunctionType):
            csum.update(strip_comments(inspect.getsource(arg)).encode())
            c = inspect.getclosurevars(arg)
            for v in {**c.nonlocals, **c.globals}.values():
                if v is not arg:
                    checksum(v, csum=csum)
#        elif hasattr(arg, '__code__'):
#            logger.debug('Checksum via __code__ for %s', str(arg))
#            bstr += arg.__code__.co_code
#            if arg.__closure__ is not None:
#                for cell in arg.__closure__:
#                    bstr += str(checksum(cell.cell_contents)).encode()
        elif isinstance(arg, functools.partial):
            logger.debug('Checksum via partial for %s', str(arg))
            checksum(arg.func, csum=csum)
            for x in arg.args:
                checksum(x, csum=csum)
            for k in sorted(arg.keywords.keys()):
                csum.update(k.encode())
                checksum(arg.keywords[k], csum=csum)
        elif isinstance(arg, np.ndarray):
            csum.update(arg.tobytes())
        else:
            logger.debug('Checksum via str for %s', str(arg))
            csum.update(str(arg).encode())

    return int.from_bytes(csum.digest(), 'big')
    