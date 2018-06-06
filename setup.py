#!/usr/bin/env python3

from distutils.core import setup
from RecBlast import __version__
setup(name='RecBlast',
      version=__version__,
      description='A flexible module to perform Reciprocal Best-Hit BLAST searches',
      author='Juan Manuel Vazquez',
      author_email='juan@vazquez.bio',
      url='https://www.vazquez.bio/RecBlast',
      install_requires=['biopython>=1.69', 'dill>=0.2.6', 'multiprocess', 'mygene'],
      packages=['RecBlast'])
