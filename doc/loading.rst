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
has to be specified with the appropiate keyword argument.
For details see :func:`mdevaluate.load_simulation`.

The function will return a coordinates object, for the whoole system.
A subset of the system may be obtained directly from the coordinates object by
calling its :func:`mdevaluate.coordinates.Coordinates.subset` method.
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

  trajectory = md.load_simulation('/data/niels/tutorial')
  CW_atoms = trajectory.subset(atom_name='CW.')

And that's it, now one can evaluate stuff for this subset of atoms.


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
This is done with the commandline tool ``index-xtc``::

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
  To this time, even though implented, the usage of trr-files has not been fully tested.

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

Note that the coordinate transformtion is not limited to center of masses at all.
Look at the definition of :func:`centers_of_mass` for hints how to implement a different transformation.
On important bit is the decorator ``@coordinates_map`` that is necessary for the transformation to work.
