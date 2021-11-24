#!/usr/bin/env python3
#  vim: set ts=4 sw=4 tw=0 :

from distutils.core import setup
from distutils.extension import Extension

import numpy

try:
    from Cython.Build import cythonize

    USE_CYTHON = True
    EXT = "pyx"
except ImportError:
    USE_CYTHON = False
    EXT = "c"
    cythonize = lambda ext: ext

extensions = [
    Extension(
        "cp2k",
        ["cp2k.{}".format(EXT)],
        include_dirs=[numpy.get_include()],
        libraries=["cp2k"],
    ),
]

setup(name="cp2k", ext_modules=cythonize(extensions))
