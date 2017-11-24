
Evaluation of static properties
===============================

.. note::
  All examples in this section assume, that the packages has been imported and a trajectory was loaded::

    import mdevaluate.distribution as dist

    coords = mdevaluate.open('/path/to/simulation')

Static properties of the system, like density distribution or pair correlation function,
can be evaluated with the :mod:`mdevaluate.distribution` module.
It provides the function :func:`mdevaluate.distribution.time_average`
that computes the average of a property over the whole trajectory.
An example call of this function is::

  tetra = dist.time_average(dist.tetrahedral_order, coords)

This will calculate the average of the tetrahedral order parameter for each atom.
The first argument of :func:`time_average` is a function that takes one argument.
It will be called for each frame in the trajectory and the output of this function
is than averaged over all these frames.

Slicing of the trajectory
-------------------------

In most cases averaging each frame of the trajectory is not necessary,
since the conformation of the atoms doesn't change significantly between two frames.
Hence it is sufficient to skip some frames without suffering significant statistics.
The exact amount of frames which can be skipped before the statistics suffer depends strongly
on the calculated property, therefore it has to be chosen manually.
For this purpose the Coordinates objects can be sliced like any python list::

  tetra = dist.time_average(dist.tetrahedral_order, coords[1000::50])

This makes it possible to skip a number of frames at the start (or end) and with every step.
The above call would start with frame 1000 of the trajectory and evaluate each 50th frame until the end.
Since the number of frames read and evaluated is reduced by about a factor of 50, the computational cost will decrease accordingly.

Calculating distributions
-------------------------

In many cases the static distributions of a property is of interest.
For example, the tetrahedral order parameter is often wanted as a distribution.
This can too be calculated with ``time_average`` but the bins of the distribution have to be specified::

  from functools import partial
  func = partial(dist.tetrahedral_order_distribution, bins=np.linspace(-3, 1, 401)
  tetra_dist = dist.time_average(func, coords)

The bins (which are ultimately used with the function :func:`numpy.histogram`) are specified
by partially evaluating the evaluation function with :func:`functools.partial`.
See the documentation of :func:`numpy.histogram` for details on bin specification.

.. note::
  If :func:`numpy.histogram` is used with :func:`time_average` the bins have to be given explicitly.
  When not specified, the bins will be chosen automatically for each call of ``histogram`` leading to
  different bins for each frame, hence an incorrect average.

Advanced evaluations
--------------------

The function that will be evaluated by ``time_average`` can return numpy arrays of arbitrary shape.
It is for example possible to calculate the distribution of a property for several subsets of the system at once::

  def subset_tetra(frame, bins):
      tetra = dist.tetrahedral_order(frame)
      return array([np.histogram(tetra[0::2], bins=bins),
                    np.histogram(tetra[1::2], bins=bins),])

  func = partial(subset, bins=np.linspace(-1,1,201))
  tetra_subdist = dist.time_average(func, coords)

In this example the tetrahedral order parameter is first calculated for each atom of the system.
Then the distribution is calculated for two subsets, containing atoms (0, 2, 4, 6, ...) and (1, 3, 5, 7, ...).
