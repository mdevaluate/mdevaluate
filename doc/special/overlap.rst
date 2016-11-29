Computing the Overlap Function
==============================

The overlap function is defined as the portion of particles of a given set,
whose positions *overlap* after a given time :math:`t` with the reference configuration at :math:`t=0`.
This is calculated as follows:
The Initial positions define spheres of a given radius :math:`r` which then are used
to test how many of the particles at a later time are found within those spheres.
Normalized by the number of spheres this gives the correlation of the configurational overlap.

.. math::

   Q(t) = \frac{1}{N} \left\langle \sum\limits_{i=1}^N n_i(t) \right\rangle

Where :math:`n_i(t)` defines the :math:`N` spheres, with :math:`n_i(t)=1` if a particle
is found within this sphere at time :math:`t` and :math:`n_i(0) = 1` for :math:`1\leq i \leq N`.

Evaluation with mdevaluate
--------------------------

Computation of the overlap requires the relatively expensive computation of next neighbor distances,
which scales with the order of :math:`\mathcal{O}(N^2)`.
There are more efficient ways for the solution of this problem, the one used here is
the so called :class:`~scipy.spatial.cKDTree`.
This is much more efficient and allows to compute the overlap relatively fast::

  OW = md.open('/path/to/sim').subset(atom_name='OW')
  tree = md.coordinates.CoordinatesKDTree(OW)
  Qol = md.correlation.shifted_correlation(
      partial(md.correlation.overlap, crds_tree=tree, radius=0.11),
      OW
  )

As seen above, mdevaluate provides the function :func:`~mdevaluate.correlation.overlap`
for this evaluation, which uses a special object of type :class:`~mdevaluate.coordinates.CoordinatesKDTree`
for the neighbor search.
The latter provides two features, necessary for the computation:
First it computes a :class:`~scipy.spatial.cKDTree` for each necessary frame of the trajectory;
second it caches those trees, since assembly of KDTrees is expensive.
The size of the cache can be controlled with the keyword argument ``maxsize`` of the CoordinatesKDTree initialization.

Note that this class uses the C version (hence the lowercase C) rather than
the pure Python version :class:`~scipy.spatial.KDTree` since the latter is significantly slower.
The only downside is, that the C version had a memory leak before SciPy 0.17,
but as long as a recent version of SciPy is used, this shouldn't be a problem.

Overlap of a Subsystem
----------------------

In many cases the overlap of a subsystem, e.g. a spatial region, should be computed.
This is done by selecting a subset of the initial configuration before defining the spheres.
The overlap is then probed with the whole system.
This has two benefits:

1. It yields the correct results
2. The KDTree structures are smaller and thereby less computation and memory expensive

An example of a spatial resolved analysis, where ``OW`` is loaded as above::

  selector = partial(
      md.coordinates.spatial_selector,
      transform=md.coordinates.spherical_radius,
      rmin=1.0,
      rmax=1.5
  )
  tree = md.coordinates.CoordinatesKDTree(OW, selector=selector)
  Qol = md.correlation.shifted_correlation(
      partial(md.correlation.overlap, crds_tree=tree, radius=0.11),
      OW
  )

This computes the overlap of OW atoms in the region :math:`1.0 \leq r \leq 1.5`.
This method can of course be used to probe the overlap of any subsystem, which is selected by the given selector function.
It should return a viable index for a (m, 3) sized NumPy array when called with original frame of size (N, 3)::

  subset = frame[selector(frame)]
