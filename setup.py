from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

include_dirs = [numpy.get_include()]

extensions = [
    Extension('mdevaluate.gromacs.coordinates', [
              'mdevaluate/gromacs/coordinates.pyx'], include_dirs=include_dirs),
    Extension('mdevaluate.gromacs.logarithmic', [
              'mdevaluate/gromacs/logarithmic.pyx'], include_dirs=include_dirs),
]


setup(
    name='mdevaluate',
    description='Collection of python utilities for md simulations',
    author_email='mbartelm@nmr.physik.tu-darmstadt.de',

    packages=['mdevaluate',
              'mdevaluate.gromacs',
              'mdevaluate.meta'],

    version='1.0.0',
    requires=['numpy', 'scipy', 'Cython'],
    ext_modules=cythonize(extensions),

    entry_points={
        'console_scripts': [
            'index-xtc=mdevaluate.gromacs.xtcindex:main',
            'demux-xtc=mdevaluate.gromacs.xtcdemux:main'
        ],
    },
)
