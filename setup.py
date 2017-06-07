#!/usr/bin/env python3

from distutils.core import setup
setup(name='recblast_MP',
      version='1.4dev0',
      description='A flexible module to perform Reciprocal Best-Hit BLAST searches',
      author='Juan Manuel Vazquez',
      author_email='juan@vazquez.bio',
      url='https://www.vazquez.bio/RecBlast',
      install_requires=['biopython>=1.69', 'dill', 'multiprocess'],
      py_modules=['recblast_MP','Auxilliary', 'test'])
