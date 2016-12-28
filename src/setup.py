# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 18:48:08 2014

@author: Weigang
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("BGM2_utils", ["BGM2_utils.pyx"],
                             include_dirs = [numpy.get_include()])]
)