# run 
# 		$ python setup_pYqRand.py build_ext --inplace
# to create shared library pYqRand.so, which can be imported into python as
# 		>>> impoprt pYqRand
# 		>>> from pYqRand import *

from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
    name = "pYqRand",
    ext_modules = cythonize('pYqRand.pyx'))
