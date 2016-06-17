Getting started
===============

The mdevaluate package has to be installed to the local python distribution.
Simply importing the modules will not work, since some of the Gromacs related
stuff has to be compiled, which is done with cython when installing the package.

The easiest way to use mdevaluate on the institutes network is to use the mdevaluate module::

  module use /data/niels/modules/modulefiles
  module load mdevaluate/1.2

Mdevaluate is provided in different versions: Numbered *release* versions (e.g. 1.2)
and the special version `dev`, which provides the latest version of the code.
The `dev` version should not be used in *production* (e.g. in the end phase of a thesis)
since it may contain broken code.

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

- Cython
- NumPy
- SciPy
- pygmx

In a future version of mdevaluate the whole gromacs dependency (and thereby the C/C++ dependencies)
will be moved to a separate package ``pygmx``.
The latest version of mdevaluate now depends partially on this new package, hence it has to be installed to.
See the `pygmx repository <https://chaos3.fkp.physik.tu-darmstadt.de/diffusion/GMX/>`_ for instructions.

Installation
------------

First obtain the source code of mdevaluate from the `git repository <https://bitbucket.org/fkp-md/mdevaluate>`_.
Navigate to the source directory and run the installation::

  python3 setup.py install

This will compile and install the package to your local python version.
