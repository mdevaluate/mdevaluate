
Tutorial
========

In this tutorial the basic usage of the python package `mdevaluate` is shown.
For further information checkout the  :doc:`modules`.
All the examples in this document should run on the local network, were ``/data`` is accessible.

Contents:

1. :ref:`tutorial-loading`
2. :ref:`tutorial-evaluation`

.. _tutorial-loading:

Loading simulation data
-----------------------

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

.. _tutorial-evaluation:

Evaluation of simulation data
-----------------------------

The evaluation of molecular dynamics siumlations can mainly be split into two parts:

1. Static properties
2. Dynamic correlations

Distributions
+++++++++++++

Static properties like radial pair distributions can be calculated with the function :func:`mdevaluate.distribution.time_average`.
This function calculates an average of the given function frame by frame over the whole trajectory.

.. Mathematically this is expressed by
.. $$
.. F(r) = \\frac{1}{N}\\sum_{i=0}^N f(R(t_i),r)
.. $$

The function that is averaged should take exactly one arguent wich is the list of coordinates in the frame.
Therefore more complex functions have to be partially evaluated, which can be done with :func:`functools.partial`::

  from mdevaluate import distribution
  from functools import partial
  from numpy import linspace

  gr = partial(distribution.radial_pair_distribution,
               bins=linspace(0,2),
               box=trajectory[0].box.diagonal())

  pair = distribution.time_average(gr, com_amim)

In the above example, the radial pair distribution between the centers of mass of the amim molecules is calculated.
The parameter ``bins`` defines the distances for which the function is computed,
the ``box`` parameter defines the periodic boundary condtions that are considered in the calculation.

Correlations
++++++++++++

Dynamic properties like mean square displacement are calculated with the
function :func:`mdevaluate.correlation.shifted_correlation`.
This function takes a correlation function and calculates the avaraged
time series of it, by shifting a time intervall over the trajectory.

::

  from mdevaluate import correlation

  ndx, msd_amim = correlation.shifted_correlation(correlation.msd, com_amim)
  time = [trajectory[i].time for i in ndx]
  msd_amim = msd_amim.mean(axis=0)
  plot(time,msd_amim)

The result of :func:`shifted_correlation` are two lists, the first one (``ndx``)
contains the indices of the frames that have been used for the correlation.
The second list ``msd_amim`` is the actual data, that is returned without avaraging over the shifted time intervalls.
As seen above the actual time steps can be read from the trajectory through the indices
and the correlation data should normally be avaraged along the first axis.
