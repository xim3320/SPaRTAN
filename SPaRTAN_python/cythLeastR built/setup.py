from distutils.core import setup, Extension
from Cython.Build import cythonize
from numpy import get_include

ext_main = Extension('cythLeastR', ['cythLeastR/cythLeastR.pyx', 'cythLeastR/ep21R.c', 'cythLeastR/eppMatrix.c', 'cythLeastR/epph.c'], include_dirs=['.', 'cythLeastR', get_include()])

setup(ext_modules=cythonize([ext_main]))
