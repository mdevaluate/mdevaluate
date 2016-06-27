Spatial Resolved Analysis
=========================

This section describes how spatially resolved correlation can be analyzed with mdevaluate.
This guide assumes that the variable ``traj`` holds a trajectory where the subset of atoms that should be analyzed are selected.
For example::

  traj = md.open('/path/to/sim', cached=1000).subset(atom_name='OW')

Which would load a simulation from the directory ``/path/to/sim`` and select all ``OW`` atoms.
Note that for this use case, the caching is quite useful since it enables us to iterate over spatial regions
without significant time penalty.
Moving on let's calculate the ISF of water oxygens with spherical radius between 0.5 and 0.7 nm::

  from functools import partial
  func = partial(md.correlation.isf, q=22.7)
  selector = partial(
    md.coordinates.spatial_selector,
    transform=md.coordinates.spherical_radius,
    rmin=0.5, rmax=0.7
  )
  t, S = md.correlation.shifted_correlation(
    func, traj,
    correlation=md.correlation.subensemble_correlation(selector)
  )

To explain how this works, let's go through the code from bottom to top.
The spatial filtering is done inside the shifted_correlation by the function
:func:`mdevaluate.correlation.subensemble_correlation`.
This function takes a selector function as argument that should take a frame as input
and return the selection of the coordinates that should be selected.
A new selection is taken for the starting frame of each shifted time segment.

In this case the selection is done with the function :func:`mdevaluate.coordinates.spatial_selector`.
This function takes four arguments, the first being the frame of coordinates which is handed by :func:`subensemble_correlation`.
The second argument is a transformation function, which transforms the input coordinates to the coordinate which will be filtered,
in this case the spherical radius.
The two last arguments define the minimum and maximum value of this quantity.
