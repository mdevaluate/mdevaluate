Spatial Resolved Analysis
=========================

This section describes how spatially resolved correlation can be analyzed with mdevaluate.
This guide assumes that the variable ``traj`` holds a trajectory where the subset of atoms that should be analyzed are selected.
For example::

  traj = md.open('/path/to/sim', cached=1000).subset(atom_name='OW')

Which would load a simulation from the directory ``/path/to/sim`` and select all ``OW`` atoms.
Note that for this use case, the caching is quite useful since it enables us to iterate over spatial regions
without significant time penalty.
Moving on let's define the function that will be correlated, in this case the ISF for water::

  from functools import partial
  func = partial(md.correlation.isf, q=22.7)

In this case, the analysis will be resolved by the spherical radius of the atoms.
We define a selector function that selects atoms between radii 0.5 and 0.7::

  selector = partial(
    md.coordinates.spatial_selector,
    radii=md.coordinates.spherical_radius(traj),
    rmin=0.5, rmax=0.7
  )

The spatial filtering is done inside the correlation, therefore we use the function
:func:`mdevaluate.correlation.subensemble_correlation`::

  t, S = md.correlation.shifted_correlation(
    func, traj,
    correlation=md.correlation.subensemble_correlation(selector)
  )
