General Hints for Python Programming
====================================

This page collects some general hints for data centered programming  with Python.
Some resources for tutorials on the topics can be found here:

* http://www.scipy-lectures.org/
* PyCon-Talk on Numpy arrays: `Losing Your Loops, by Jake VanderPlas <https://www.youtube.com/watch?v=EEUXKG97YRw>`_

Programming Environments
------------------------

There exist different environments for Python programming, each with their pros and cons.
Some examples are:

* **IPython Console**: The most basic way to use Python is on the interactive console, the ipython console is a suffisticated Python console. After the mdevaluate module is loaded, ipython can be started with the command ``ipython``.
* **Jupyter Notebook**: Provides a Mathematica-style notebook, which is accesed through a web browser. After the mdevaluate module is loaded a (local) notebook server can be started with the command ``jupyter-notebook``. See the help menu in the notebook for a short introduction and http://jupyter.org/ for a detailed user guide.
* **Atom Editor**: When developing more complex code, like modules an editor comes in handy. Besides basic preinstalled editors (e.g. Gedit) the `atom editor <https://atom.io>`_ is a nice option. Recommended atpm packages for Python development are: language-python, autocomplete-python and linter-flake8.

Common Pitfalls
---------------

* **For-Loops**: The biggest pitfall of data-intensive Python programming are ``for``-loops. Those loops perform bad in Python, but can be avoided in most cases through Numpy arrays, see the mentioned talk by Jake VdP.
* **Non-Portable Code**: Most non-programmers tend to write complex scripts. It's always advisable to source out your code into seperate Python modules (i.e. seperate files) and split the code into reusable functions. Since these modules can be imported from any Python code, this will save time in the long run and often reduces errors.


Pandas Dataframes
-----------------

Most data in Mdevaluate is handled as Numpy arrays.
For example the function :func:`~mdevaluate.correlation.shifted_correlation` returns a multidimensional array, which contains the time steps and the value of the correlation function.
As pointed out above, those arrays a good for computation and can be used to plot data with, e.g. matplotlib.
But often there is metadata associated with this data, for example the temperature or the specific subset of atoms that were analyzed.
This is the point where **`Pandas dataframes <http://pandas.pydata.org/>`_** come in handy.
Dataframes are most basically tables of samples, with named columns.
The dataframe class allows easy acces of columns by label and complex operations, like grouping by columns or merging different datasets.

As an example say we have simulations at some temperatures and want to calculate the ISF and do a KWW-Fit for each of these trajectories.
Details of the analysis will be explained at a later point of this document, thereby they will be omitted here::

  import pandas
  datasets = []

  for T in [250, 260, 270, 280, 290, 300]:
      # calculate the isf for this temperature
      t, Sqt = ...

      # DataFrames can be created from dictionaries
      datasets.append(pandas.DataFrame({'time': t, 'Sqt': Sqt, 'T': T}))

  # join the individual dataframes into one
  isf_data = pandas.concat(datasets)

  # Now calculate the KWW fits for each temperature
  from scipy.optimize import curve_fit
  from mdevaluate.functions import kww
  kww_datasets = []
  # The isf data is grouped by temperature,
  # that is the loop iterates over all T values and the part of the data where isf_data['T'] == T
  for T, data in isf_data.groupby('T'):
      fit, cuv = curve_fit(kww, data['time'], data['Sqt'])
      # DataFrames can also be cerated from arrays and a defintion of columns
      df = pandas.DataFrame(fit, columns=['A', 'τ', 'β'])
      # columns can be added dynamically
      df['T'] = T
      kww_datasets.append(df)
  kww_data = pandas.concat(kww_datasets)

  # We have two dataframes now, one with time series of the ISF at each temperature
  # and one with the fit parameters of the KWW for each temperature

  # We can merge this data into one dataframe, by the overlapping columns (i.e. 'T' in this example)
  data = pandas.merge(isf_data, kww_data)
  # We can now compute the kww fit value of each sample point of the isf in one line:
  data['kww_fit'] = kww(data['time'], data['A'], data['τ'], data['β'])
  # And plot the data, resolved by temperature.
  for T, df in data.groupby('T'):
      plot(df['time'], df['Sqt'], 'o') # The actual correlation value
      plot(df['time'], df['kww_fit'], '-') # The kww fit
