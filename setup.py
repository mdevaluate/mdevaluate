from setuptools import setup
from mdevaluate import __version__

setup(
    name='mdevaluate',
    description='Collection of python utilities for md simulations',
    author_email='niels.mueller@physik.tu-darmstadt.de',

    packages=['mdevaluate',
              'mdevaluate.meta'],

    version=__version__,
    requires=['numpy', 'scipy'],
    zip_safe=False
)
