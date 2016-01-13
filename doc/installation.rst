Installation
============

The mdevaluate package has to be installed to the local python distribution.
Simply importing the modules will only work partly, since some of the gromacs related
stuff has to be compiled, which is done with cython when installing the package.

The easiest way to use mdevaluate on the institutes network is to use an python installation from the network.
One python installation, with a recent version of mdevaluate installed is found here::

  /autohome/niels/miniconda3

To use this python installation the python (or ipython) executable has to be used.
The simplest way to achieve that is to include it in the ``PATH`` variable,
which is normally done ba adding the follwoing line to the file ``~/.bashrc``::

  export PATH=/autohome/niels/miniconda3/bin:$PATH

After doing so, the python executables will point to this directory and python or
IPython (Jupyter) may be used as usual.
The above described python installation provides Python 3 and all relevant packages
for scientific computing.
If any packages are missing, contact `Niels Müller <mailto:niels@nmr.physik.tu-darmstadt.de>`_.

Manual Installation
+++++++++++++++++++

The second option is to install the package manually to a local python installation.
Before installing mdevaluate, some python packges need to be installed on the system.
It should also be emphasized that mdevaluate will only run in **Python 3**.

Requirements
------------

The package depends on several python packages that can all be installed via pip or conda:

- Python 3
- Cython
- NumPy
- SciPy

Additional requirements for ``tudplot``:

- Matplotlib
- Seaborn

Installation
------------

First obtain the source code of mdevaluate from the `git repository <https://bitbucket.org/fkp-md/mdevaluate>`_.
Navigate to the source directory and run the installation::

  python3 setup.py install

This will compile and install the package to your local python version.
