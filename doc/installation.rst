Installation
============

The mdevaluate package has to be installed to the local python distribution.
Simply importing the modules will not work, since some of the gromacs related
stuff has to be compiled, which is done with cython when installing the package.

The easiest way to use mdevaluate on the institutes network is to use the anaconda3 module::

  module load anaconda3

This will provide a recent python 3.5 environment.
Addtionally the mdevaluate package has to be added to the `PYTHONPATH`, which may be done in your `.bashrc` file::

  export PYTHONPATH=/autohome/niels/.local/lib/python3.5/site-packages:$PYTHONPATH

The above described python installation provides Python 3 and all relevant packages
for scientific computing.
If any package is missing, you can install it locally with::

  pip install PACKAGE --user


Manual Installation
+++++++++++++++++++

The second option is to install the package manually to a local python installation.
Before installing mdevaluate, some python packges need to be installed on the system.
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
The latest version of mdevaluate now depends partailly on this new package, hence it has to be installed to.
See the `pygmx repository <https://bitbucket.org/fkp-md/pygmx>`_ for instructions.

Installation
------------

First obtain the source code of mdevaluate from the `git repository <https://bitbucket.org/fkp-md/mdevaluate>`_.
Navigate to the source directory and run the installation::

  python3 setup.py install

This will compile and install the package to your local python version.
