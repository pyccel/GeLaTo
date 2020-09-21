# -*- coding: UTF-8 -*-
#! /usr/bin/python

from pathlib import Path
from setuptools import setup, find_packages

# ...
# Read library version into '__version__' variable
path = Path(__file__).parent / 'gelato' / 'version.py'
exec(path.read_text())
# ...

NAME    = 'gelato'
VERSION = __version__
AUTHOR  = 'Ahmed Ratnani'
EMAIL   = 'ahmed.ratnani@ipp.mpg.de'
URL     = 'https://github.com/ratnania/GeLaTo'
DESCR   = 'A Python library for Generalized Locally Toeplitz theory for Isogeometric Analysis.'
KEYWORDS = ['math']
LICENSE = "LICENSE"

setup_args = dict(
    name                 = NAME,
    version              = VERSION,
    description          = DESCR,
    long_description     = open('README.rst').read(),
    author               = AUTHOR,
    author_email         = EMAIL,
    license              = LICENSE,
    keywords             = KEYWORDS,
    url                  = URL,
#    download_url     = URL+'/tarball/master',
)

install_requires = [

    # Third-party libraries from PyPi
    'sympy>=1.2',
    'numpy>=1.13',
    'scipy>=0.18',
    'matplotlib',

    # Our libraries from PyPi
    'sympde>=0.10',
]

# ...
packages = find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"])
# ...

def setup_package():

    setup(packages = packages,
          include_package_data = True,
          install_requires = install_requires,
          zip_safe = True,
          **setup_args)


if __name__ == "__main__":
    setup_package()
