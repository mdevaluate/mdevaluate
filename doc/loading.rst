Loading of simulation data
==========================

.. note::
  The process of loading simulations has been simplified in the latest version of mdevaluate.
  The previous way is described below for completeness and people with an older code base.

Mdevaulate provides a convenient function :func:`mdevaluate.load_simulation`
that loads a simulation more or less automatically.
It takes a path as input and looks for all files it needs in this directory.

For information about the topology either a `tpr` or `gro` a file is read,
where the former is the preferred choice.
Trajectory data will be read from a xtc file.
If the directory contains more than one file of any type, the desired file
has to be specified with the appropriate keyword argument.
For details see :func:`mdevaluate.load_simulation`.

The function will return a coordinates object, for the whole system.
A subset of the system may be obtained directly from the coordinates object by
calling its :func:`~mdevaluate.coordinates.Coordinates.subset` method.
This function accepts the same input as :func:`mdevaluate.atoms.AtomSubset.subset`.
A new feature that was introduced in the function is the possibility to chose
atoms with regular expressions.

Example
-------

The following code loads the example trajectory and selects a subset of all CW atoms.
Since there are two CW atoms in each molecule (CW1 and CW2) a regular expression is
used when selecting the subset.

::

  import mdevaluate as md

  trajectory = md.open('/data/niels/tutorial')
  CW_atoms = trajectory.subset(atom_name='CW.')

And that's it, now one can evaluate stuff for this subset of atoms.

Selecting a subset
------------------

As shown in the example above it is often necessary to select a subset of the system for analysis.
This can be a special group of atoms (e.g. all C atoms) or a whole residue for which the center of mass should be computed.
Subsets are selected with the :func:`~mdevaluate.Coordinates.subset` method of Coordinates objects.

This method accepts four keyword arguments, with which the atom name, residue name and residue id or atom indices can be specified.
The former two name arguments accept a regular expression which allows two include several different names in one subset.
Some examples:

- All carbon atoms (which are named CW1, CT1, CA, ...): ``tr.subset(atom_name='C.*')``
- Atoms NA1, NA2 and OW: ``tr.subset(atom_name='NA.|OW')``
- All oxygen atoms of residue EG: ``tr.subset(atom_name='O.*', residue_name='EG')``


Specifying data files
---------------------

The above example only works if the directory contains exactly one tpr file and
one xtc file.
If your data files are located in subdirectories or multiple files of these types exist,
they can be specified by the keywords ``topology`` and ``trajectory``.
Those filenames can be a relative path to the simulation directory and can also make
use of *shell globing*. For example::

  traj = md.open('/path/to/sim', topology='atoms.gro', trajectory='out/traj_*.xtc')

Note that the topology can be specified as a gro file, with the limitation that
only atom and residue names will be read from those files.
Information about atom masses and charges for example will only be read from tpr files,
therefore it is generally recommended to use the latter topologies.

The trajectory above is specified through a shell globing, meaning the ``*`` may be
expanded to any string (without containing a forward slash).
If more than one file exists which match this pattern an error will be raised,
since the trajectory can not be identified clearly.

Caching of frames
-----------------

One bottleneck in the analysis of MD data is the reading speed of the trajectory.
In many cases frames will be needed repeatedly and hence the amount of time spend reading
data from disk (or worse over the network) is huge.
Therefore the mdevaluate package implements a simple caching mechanism, which holds
on to a number of read frames.
The downside if this is increased memory usage which may slow down the computation too.

Caching is done on the level of the trajectory readers, so that all ``Coordinate`` and
``CoordinateMap`` objects working on the same trajectory will be sharing a cache.
Caching has to be activated when opening a trajectory::

  traj = md.open('/path/to/sim', cached=True)

The ``cached`` keyword takes either a boolean, a integer or None as input value.
The value of ``cached`` controls the size of the cache and thereby the additional memory usage.
Setting it to True will activate the caching with a maximum size of 128 frames,
with an integer any other maximum size may be set.
The special value ``None`` will set the cache size to infinite, so all frames will be cached.
This will prevent the frames from being loaded twice but can also consume a whole lot of memory,
since a single frame can easily take 1 MB of memory.

Freeing cached frames
+++++++++++++++++++++

In some scenarios it may be advisable to free cached frames which are no longer needed.
For this case the reader has a function ``clear_cache()``.
The current state of the cache can be displayed with the ``cache_info`` property::

  >>> traj.frames.cache_info
  CacheInfo(hits=12, misses=20, maxsize=128, currsize=20)
  >>> traj.frames.clear_cache()
  >>> traj.frames.cache_info
  CacheInfo(hits=0, misses=0, maxsize=128, currsize=0)

Clearing the cache when it is not needed anymore is advisable since this will help the
Python interpreter to reuse the memory.


The deprecated way of loading simulations
-----------------------------------------

For the evaluation, two data sets are necessary:

1. Information about the atoms of the simulation, which is read from gro-files
2. The trajectory of the simulation, which is read from xtc- or trr-files

These two data sets are combined to a coordinates object, that selects only the requested atoms.
This object can also be used to perform coordinate transformations like computing a center of mass.

Atoms
+++++

Atom information is read from gro-files via the function :func:`mdevaluate.atoms.from_grofile`,
which takes the path to a gro-file as argument and an optional gromacs index file.

::

  from mdevaluate import atoms

  all_atoms = atoms.from_grofile('/data/niels/tutorial/conf.gro',
                                 index='/data/niels/tutorial/index.ndx')

For many evaluations, a subset of the system must be selected.
Atoms can be selected by name, residue or direct indices::

  H11_atoms = all_atoms.subset(atom_name='H11')
  amim_atoms = all_atoms.subset(residue_name='AMIM')

Atom subset can be combined with logical operators to obtain the intersection or union of two or more subsets.
The union is equal to a logical or, the intersection is equal to a logical and::

  union = H11_atoms | amim_atoms
  intersection = H11_atoms & amim_atoms

Excluding atoms from a larger subset can also be done easily with negation::

  exclusion = amim_atoms & ~H11_atoms

Trajectory
++++++++++

The trajectory is read from xtc-files or trr-files, usually the former are used.
For effective data loading, these files have to be indexed **once** before they can be used with mdevaluate.
This is done with the command line tool ``index-xtc``::

  $ index-xtc /data/niels/tutorial/traj.xtc

The trajectory is than read with an appropiate reader::

  from mdevaluate.gromacs.reader import XTCReader

  trajectory = XTCReader('/data/niels/tutorial/traj.xtc')

From this trajectory, single frames can be selected by index::

  frame = trajectory[42]
  print(frame.time)
  print(frame.box)
  print(frame.coordinates)

.. warning::
  To this time, even though implemented, the usage of trr-files has not been fully tested.

Coordinates
+++++++++++

A coordinates object for the evaluation is obtained by combining the trajectory and an atom subset::

  from mdevaluate import coordinates

  cords_amim = coordinates.Coordinates(trajectory, atom_subset=amim_atoms)

These coordinates can be transformed if necessary.
The center of mass can be computed with the function :func:`mdevaluate.coordinates.centers_of_mass`,
which takes coordinates and a list of masses as input.
Only the first set of atom masses has to be given, which will be repeated for the rest of the atoms.
The only requirement is that the length of the list of atoms is an integral multiple of the number of masses given.

To compute a center of mass of the amim molecule::

  masses = [14]*2 + [1]*15 + [12]*8
  com_amim = coordinates.centers_of_mass(cords_amim, masses=masses)

Note that the coordinate transformation is not limited to center of masses at all.
Look at the definition of :func:`centers_of_mass` for hints how to implement a different transformation.
On important bit is the decorator ``@coordinates_map`` that is necessary for the transformation to work.
