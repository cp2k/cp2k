#!/usr/bin/env python
import sys
if(len(sys.argv)==2 and sys.argv[-1]=="--selftest"):
       sys.exit(0)

from setuptools import setup

setup(name='fprettify',
      description='auto-formatter for modern fortran source code',
      author='Mohamed Fawzi, Patrick Seewald, Ole Schuett',
      license = 'GPL',
      entry_points={'console_scripts': ['fprettify = fprettify:main']},
      py_modules=['fprettify',
          'formatting.normalizeFortranFile',
          'formatting.reformatFortranFile',
          'formatting.replacer',
          'formatting.selftest']
     )
