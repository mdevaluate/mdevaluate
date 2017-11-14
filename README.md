# mdevaluate

Mdevaluate is a Python package to perform analyses of  Molecular Dynamics simulations.
An online documentation is available at [mdevaluate.github.io](https://mdevaluate.github.io).
Mdevaluate provides a flexible interface for the detailed analysis of dynamical and statical properties of molecular systems.
It's main focus is the analysis of Gromacs data, but with the help of external packages ([MDAnalysis](https://www.mdanalysis.org/)) 
it can also handle file formats, used by other simulation software.
	
	import mdevaluate as md

	# load the simulation
	tr = md.open(
		directory='/path/to/simulation',
		topology='topol.tpr',
		trajectory='traj.xtc'
	)
	# select a subset of atoms
	water_oxygen = tr.subset(residue_name='SOL', atom_name='OW')

	# calculate the mean squared displacement for this subset
	time, msd = md.correlation.shifted_correlation(
		md.correlation.msd, water_oxygen, average=True
	)


## Installation

The package requires the Python package [pygmx](https://github.com/mdevaluate/pygmx),
which handles reading of Gromacs file formats.
Installation of pygmx is described in its own repository.

The mdevaluate package itself is plain Python code and, hence, can be imported from its directory directly, 
or may be installed via setuptools to the local Python environment by running
	
	python setup.py install


## Running the tests

Mdevaluate includes a test suite that can be used to check if the installation was succesful.
It is based on `py.test` and located in the test directory of this repository.
Make sure py.test is installed and run `py.test` within the repository to check if all tests pass.


