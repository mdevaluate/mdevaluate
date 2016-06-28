Automatic Saving of Analysis Data
=================================

Mdevaluate provides a functionality to save the result of analysis functions automatically.
The data is saved to a file after it was computed.
If an analysis was done in the exact same way before, the result is loaded from this file.

This function may be activated through calling :func:`mdevaluate.autosave.enable`, which takes a directory as input.
If this directory is a relative path (e.g. no trailing slash) the results will be saved in a location
relative to the directory of the trajectory file.
If the output files of your simulations are located in a subdirectory, like ``/path/to/sim/Output`` it is possible
to specify the auto save location as ``../data`` such that the result files will be placed under ``/path/to/sim/data``.

At the moment the two functions which use this behavior are:

- :func:`~mdevaluate.correlation.shifted_correlation`
- :func:`~mdevaluate.distribution.time_average`

Any other function can make use of the autosave mechanism by decorating it with :func:`mdevaluate.autosave.autosave_data`.

A full example
--------------

This is how it works, for a detailed explanation see below::

  import mdevaluate as md
  md.autosave.enable('data')
  water = md.open('/path/to/sim').subset(atom_name='OW')
  md.correlation.shifted_correlation(
      md.correlation.msd,
      water,
      description='test'
  )
  # The result will be saved to the file:
  # /path/to/sim/data/shifted_correlation_msd_OW_test.npz

Checksum of the Analysis Call
-----------------------------

The autosave module calculates a checksum for each call of an analysis function,
which is used to validate a present the data file.
This way the result should only be loaded from file if the analysis is exactly the same.
This includes the function code that is evaluated, so the result will be recomputed if any bit of the code changes.
But there is always the possibility that checksums coincide accidentally,
by chance or due to a bug in the code, which should be kept in mind when using this functionality.

Special Keyword Arguments
-------------------------

The autosave module introduces two special keyword arguments to the decorated functions:

- ``autoload``: This prevents the loading of previously calculated data even if a valid file was found.
- ``description``: A descriptive string of the specific analysis, see below.

Those keywords may be passed to those function (shifted_correlation, time_average) like any other keyword argument.
If autosave was not enabled, they will be ignored.

File names and Analysis Descriptions
------------------------------------

The evaluated data is saved to human readable files, whose name is derived from the function call
and the automatic description of the subset.
The latter one is assigned based on the ``atom_name`` and ``residue_name`` of the :func:`~mdevaluate.atoms.AtomSubset.subset` method.

In some cases this is not enough, for example if the same subset is analyzed spatially resolved,
which would lead to identical filenames that would be overwritten.
Therefore a more detailed description of each specific analysis call needs to be provided.
For this reason the autosave module introduces the before mentioned keyword argument ``description``.
The value of this keyword is appended to the filename and in addition if any of
the other arguments of the function call has a attribute description, this will appended as well.
For example this (pseudo) code will lead to the filename ``shifted_correlation_isf_OW_1-2nm_nice.npz``::

  OW = traj.subset(atom_name='OW')

  corr = subensemble_correlation(spatial_selector)
  corr.description = '1-2nm'

  shifted_correlation(
    isf,
    OW,
    correlation=corr,
    description='nice'
  )


Reusing the autosaved data
--------------------------

The results of the functions are saved in NumPy's  npz format, see :func:`numpy.savez`.
If the result should be used in a different place, it can either be loaded with
:func:`numpy.load` or :func:`mdevaluate.autosave.load_data`.
The latter function will return the result of the function call directly, the former
returns a dict with the keys ``checksum`` and ``data``, the latter yielding the results data.
