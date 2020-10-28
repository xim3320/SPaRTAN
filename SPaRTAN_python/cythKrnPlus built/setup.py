# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 19:44:19 2020

@author: XIM33
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(name='kronecker',
      ext_modules=cythonize("cythkrnPlus.pyx"))
	  
