from setuptools import setup
import os.path

users = [
    'inken',
    'kai',
    'kbeck',
    'kurt',
    'lisa',
    'marius',
    'mbartelm',
    'niels',
    'tamisra',
    'wieth',
]


packages = [
    'mdevaluate_user.' + user for user in users
]


setup(
    name='mdevaluate_user',
    description='Collection of python utilities for md simulations',
    author_email='mbartelm@nmr.physik.tu-darmstadt.de',
    packages = packages,
    version='1.0.0',
)