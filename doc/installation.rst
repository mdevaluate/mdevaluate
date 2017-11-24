Installation
============

Mdevaluate itself is a pure Python package and can be imported directly from the source directory, if needed.
The Gromacs dependency pygmx has to be installed into the Python distribution,
since parts are compiled with Cython.

Requirements
------------

The package depends on some python packages that can all be installed via pip or conda:

- Python 3.5 (or higher)
- NumPy
- SciPy


Install pygmx & mdevaluate
---------------------------

To instal pygmx, first get the source from its repository, https://github.com/mdevaluate/pygmx.
Installation instructions are given in the respective readme file.
Two steps have to be performed:

1. Install Gromacs 2016
2. Install pygmx 

When this requirement is met, installing mdevaluate simply means getting the source code from the repository and running 
	
	python setup.py install 

form within the source directory.


Running Tests
-------------

Some tests are included with the source that can be used too test the installation.
The testsuite requires `pytest <https://pytest.org>`_.
To run the test simply execute

	pytest

in the source directory.