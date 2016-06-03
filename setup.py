from setuptools import setup

setup(
    name='mdevaluate',
    description='Collection of python utilities for md simulations',
    author_email='niels.mueller@physik.tu-darmstadt.de',

    packages=['mdevaluate',
              'mdevaluate.meta'],

    version='1.2',
    requires=['numpy', 'scipy'],
)
