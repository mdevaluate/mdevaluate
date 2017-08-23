import os
import numpy as np
import functools

from .checksum import checksum
from .logging import logger

autosave_directory = None
load_autosave_data = False
verbose_print = True
user_autosave_directory = os.path.join(os.environ['HOME'], '.mdevaluate/autosave')


def notify(msg):
    if verbose_print:
        logger.info(msg)
    else:
        logger.debug(msg)


def enable(dir, load_data=True, verbose=True):
    """
    Enable auto saving results of functions decorated with :func:`autosave_data`.

    Args:
        dir: Directory where the data should be saved.
        load_data (opt., bool): If data should also be loaded.
    """
    global autosave_directory, load_autosave_data, verbose_print
    verbose_print = verbose
    # absolute = os.path.abspath(dir)
    # os.makedirs(absolute, exist_ok=True)
    autosave_directory = dir
    load_autosave_data = load_data
    notify('Enabled autosave in directory: {}'.format(autosave_directory))


def disable():
    """Disable autosave."""
    global autosave_directory, load_autosave_data
    autosave_directory = None
    load_autosave_data = False


class disabled:
    """
    A context manager that disbales the autosave module within its context.

    Example:
        import mdevaluate as md
        md.autosave.enable('data')
        with md.autosave.disabled():
            # Autosave functionality is disabled within this context.
            md.correlation.shifted_correlation(
                ...
            )

        # After the context is exited, autosave will work as before.
    """
    def __enter__(self):
        self._autosave_directory = autosave_directory
        disable()

    def __exit__(self, *args):
        enable(self._autosave_directory)


def get_directory(reader):
    """Get the autosave directory for a trajectory reader."""
    outdir = os.path.dirname(reader.filename)
    savedir = os.path.join(outdir, autosave_directory)
    if not os.path.exists(savedir):
        try:
            os.makedirs(savedir)
        except PermissionError:
            pass
    if not os.access(savedir, os.W_OK):
        savedir = os.path.join(user_autosave_directory, savedir.lstrip('/'))
        logger.info('Switched autosave directory to {}, since original location is not writeable.'.format(savedir))
    os.makedirs(savedir, exist_ok=True)
    return savedir


def get_filename(function, checksum, description, *args):
    """Get the autosave filename for a specific function call."""
    func_desc = function.__name__
    for arg in args:
        if hasattr(arg, '__name__'):
            func_desc += '_{}'.format(arg.__name__)
        elif isinstance(arg, functools.partial):
            func_desc += '_{}'.format(arg.func.__name__)

        if hasattr(arg, 'frames'):
            savedir = get_directory(arg.frames)

        if hasattr(arg, 'description') and arg.description != '':
            description += '_{}'.format(arg.description)
    filename = '{}_{}.npz'.format(func_desc.strip('_'), description.strip('_'))
    return os.path.join(savedir, filename)


def verify_file(filename, checksum):
    """Verify if the file matches the function call."""
    file_checksum = 0
    if os.path.exists(filename):
        data = np.load(filename)
        if 'checksum' in data:
            file_checksum = data['checksum']
    return file_checksum == checksum


def save_data(filename, checksum, data):
    """Save data and checksum to a file."""
    notify('Saving result to file: {}'.format(filename))
    try:
        data = np.array(data)
    except ValueError:
        arr = np.empty((len(data),), dtype=object)
        arr[:] = data
        data = arr

    np.savez(filename, checksum=checksum, data=data)


def load_data(filename):
    """Load data from a npz file."""
    notify('Loading result from file: {}'.format(filename))
    fdata = np.load(filename)
    if 'data' in fdata:
        return fdata['data']
    else:
        data = tuple(fdata[k] for k in sorted(fdata) if ('arr' in k))
        save_data(filename, fdata['checksum'], data)
        return data


def autosave_data(nargs, kwargs_keys=None, version=None):
    """
    Enable autosaving of results for a function.

    Args:
        nargs: Number of args which are relevant for the calculation.
        kwargs_keys (opt.): List of keyword arguments which are relevant for the calculation.
        version (opt.): 
            An optional version number of the decorated function, which replaces the checksum of
            the function code, hence the checksum does not depend on the function code.
    """
    def decorator_function(function):
        @functools.wraps(function)
        def autosave(*args, **kwargs):
            description = kwargs.pop('description', '')
            autoload = kwargs.pop('autoload', True) and load_autosave_data
            if autosave_directory is not None:
                relevant_args = list(args[:nargs])
                if kwargs_keys is not None:
                    for key in kwargs_keys:
                        if key in kwargs:
                            relevant_args.append(kwargs[key])

                if version is None:
                    csum = legacy_csum = checksum(function, *relevant_args)
                else:
                    csum = checksum(version, *relevant_args)
                    legacy_csum = checksum(function, *relevant_args)

                filename = get_filename(function, csum, description, *relevant_args)
                if autoload and (verify_file(filename, csum) or verify_file(filename, legacy_csum)):
                    result = load_data(filename)
                else:
                    result = function(*args, **kwargs)
                    save_data(filename, csum, result)
            else:
                result = function(*args, **kwargs)

            return result
        return autosave
    return decorator_function
