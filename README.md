# mdevaluate

Mdevaluate provides tools for trajectory analysis of MD simulations.
[A html documentation is available on-line](http://element.fkp.physik.tu-darmstadt.de/~niels/mdevaluate).
This documentation can also be build locally with [spinx](http://www.sphinx-doc.org),
by running `make html` in the `doc` directory.
Please refer to this documentation for a detailed description of the usage of this package.

## Setup

The package requires the Python package [pygmx](https://chaos3.fkp.physik.tu-darmstadt.de/diffusion/GMX/),
which handles reading of Gromacs file formats.
Installation of pygmx is described in its repository.

When pygmx is properly installed clone this repository and run `python3 setup.py install` to install mdevaluate.
The package can also be imported right out of the repository, by adding it to the `PYTHONPATH` variable.
