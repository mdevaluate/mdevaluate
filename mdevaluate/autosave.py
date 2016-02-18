import os
import types
import numpy as np
from .utils import merge_hashes, hash_anything as _hash
import functools
autosave_directory = None
load_autosave_data = False


def enable(dir, load_data=True):
    """
    Enable auto saving results of functions decorated with @autosave_data.

    Args:
        dir: Directory where the data should be saved.
        load_data (opt., bool): If data should also be loaded.
    """
    global autosave_directory, load_autosave_data
    absolute = os.path.abspath(dir)
    os.makedirs(absolute, exist_ok=True)
    autosave_directory = absolute
    load_autosave_data = load_data


def get_filename(function, checksum, description):
    """Get the autosave filename for a specific function call."""
    filename = '{}_{}.{}.npy'.format(function.__name__, description, checksum)
    return os.path.join(autosave_directory, filename)


def checksum(function, *args):
    """Get the checksum of a function call."""
    hashes = [_hash(function.__code__)]
    for arg in args:
        if isinstance(arg, types.FunctionType):
            hashes.append(_hash(arg.__code__))
        elif isinstance(arg, functools.partial):
            hashes.append(checksum(arg.func, *arg.args, *arg.keywords.values()))
        else:
            hashes.append(_hash(arg))
    return merge_hashes(*hashes)


def verify_file(filename, checksum):
    """Verify if the file matches the function call."""
    return os.path.exists(filename)


def autosave_data(nargs, kwargs_keys=None):
    """
    Enable autosaving of results for a function.

    Args:
        nargs: Number of args which are relevant for the calculation.
        kwargs (opt.): List of keyword arguments which are relevant for the calculation.
    """
    def decorator_function(function):
        @functools.wraps(function)
        def autosave(*args, **kwargs):
            if autosave_directory is not None:
                relevant_args = args[:nargs]
                description = kwargs.pop('description', '')
                if kwargs_keys is not None:
                    relevant_args += [kwargs[key] for key in kwargs_keys]

                csum = checksum(function, *relevant_args)
                filename = get_filename(function, csum, description)
                if load_autosave_data and verify_file(filename, csum):
                    print('Loading result from file: {}'.format(filename))
                    result = np.load(filename)
                else:
                    result = function(*args, **kwargs)
                    print('Saving result to file: {}'.format(filename))
                    np.save(filename, result)
            else:
                result = function(*args, **kwargs)

            return result
        return autosave
    return decorator_function
