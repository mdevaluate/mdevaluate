from setuptools import setup


def get_version(module):
    version = ''
    with open(module) as f:
        for line in f:
            if '__version__' in line:
                version = line.split('=')[-1].strip("' \n\t")
                break
    return version


setup(
    name='mdevaluate',
    description='Collection of python utilities for md simulations',
    author_email='niels.mueller@physik.tu-darmstadt.de',

    packages=['mdevaluate',],
    entry_points={
        'console_scripts': [
            'index-xtc = mdevaluate.cli:run'
        ]
    },
    version='dev',
    requires=['numpy', 'scipy'],
)
