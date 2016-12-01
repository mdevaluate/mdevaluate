from setuptools import setup


def get_version(module):
    version = ''
    with open(module) as f:
        for line in f:
            if '__version__' in line:
                version = line.split('=')[-1].strip("' \n\t")
                break
    return version

__version__ = get_version('mdevaluate/__init__.py')

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
