Evaluation of dynamic properties
================================

Dynamic properties like mean square displacement are calculated with the
function :func:`mdevaluate.correlation.shifted_correlation`.
This function takes a correlation function and calculates the averaged
time series of it, by shifting a time interval over the trajectory.

::

  from mdevaluate import correlation

  time, msd_amim = correlation.shifted_correlation(correlation.msd, com_amim, average=True)
  plot(time,msd_amim)

The result of :func:`shifted_correlation` are two lists, the first one (``time``)
contains the times of the frames that have been used for the correlation.
The second list ``msd_amim`` is the correlation function at these times.
If the keyword ``average=False`` is given, the correlation function for each shifted
time window will be returned.

Arguments of ``shifted_correlation``
------------------------------------

The function :func:`mdevaluate.correlation.shifted_correlation` accepts several keyword arguments.
With those arguments, the calculation of the correlation function may be controlled in detail.
The mathematical expression for a correlation function is the following:

.. math:: S(t) = \frac{1}{N} \sum_{i=1}^N C(f, R, t_i, t)

Here :math:`S(t)` denotes the correlation function at time t, :math:`R` are the coordinates of all atoms
and :math:`t_i` are the onset times (:math:`N` is the number of onset times or time windows).
Note that the outer sum and division by :math:`N` is only carried out if ``average=True``.
The onset times are defined by the keywords ``segments`` and ``window``, with
:math:`N = segments` and :math:`t_i = \frac{ (1 - window) \cdot t_{max}}{N} (i - 1)` with the total simulation time :math:`t_{max}`.
As can be seen ``segments`` gives the number of onset times and ``window`` defines the part of the simulation time the correlation is calculated for,
hence ``window - 1`` is the part of the simulation the onset times a distributed over.


:math:`C(f, R, t_0, t)` is the function that actually correlates the function :math:`f`.
For standard correlations the functions :math:`C(...)` and :math:`f` are defined as:

.. math:: C(f, R, t_0, t) = f(R(t_0), R(t_0 + t))

.. math:: f(r_0, r) = \langle s(r_0, r) \rangle

Here the brackets denote an ensemble average, small :math:`r` are coordinates of one frame and :math:`s(r_0, r)` is the value that is correlated,
e.g. for the MSD :math:`s(r_0, r) = (r - r_0)^2`.

The function :math:`C(f, R, t_0, t)` is specified by the keyword ``correlation``, the function :math:`f(r_0, r)` is given by ``function``.
