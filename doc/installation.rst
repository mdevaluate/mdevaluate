Getting started
===============

The mdevaluate package has to be installed to the local python distribution.
Simply importing the modules will not work, since some of the Gromacs related
stuff has to be compiled, which is done with cython when installing the package.

The easiest way to use mdevaluate on the institutes network is to use the mdevaluate module::

  module load mdevaluate/16.09

The mdevaluate module is provided in different versions.
Since there is no real *release schedule*, the version numbers are arbitrary and
follow a `Calendar Versioning scheme <http://calver.org>`_.
The idea of these versions is to provide a fixed version of the code from time to time,
which will not be changed, except maybe for major bugs, and assure that scripts don't stop working.
The special version `dev` provides the latest version of the code.
It should not be used in *production* (e.g. in the end phase of a thesis)
since it may contain broken code or function signatures may change.

The module will provide a recent python 3.5 environment through the anaconda3 module
with all relevant packages for scientific computing already installed.
After loading th module you can for example start a notebook by running::

  jupyter-notebook

If any package is missing, it can be installed locally with::

  pip install PACKAGE --user

Manual Installation
+++++++++++++++++++

The second option is to install the package manually to a local python installation.
Before installing mdevaluate, some python packages need to be installed on the system.
It should also be emphasized that mdevaluate will only run in **Python 3.5**.

Requirements
------------

The package depends on several python packages that can all be installed via pip or conda:

- NumPy
- SciPy

The most important dependency is ``pygmx``, which uses the Gromacs C-library to read data files.
See the `pygmx repository <https://chaos3.fkp.physik.tu-darmstadt.de/diffusion/GMX/>`_ for installation instructions.

Installation
------------

First obtain the source code of mdevaluate from the `git repository <https://chaos3.fkp.physik.tu-darmstadt.de/diffusion/MDE/>`_.
Navigate to the source directory and run the installation::

  python3 setup.py install

This will compile and install the package to your local python distribution.
