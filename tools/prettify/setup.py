#!/usr/bin/env python

from setuptools import setup

setup(name='fprettify',
      description='auto-formatter for modern fortran source code',
      author='Mohamed Fawzi, Patrick Seewald',
      license = "GPL",
      entry_points={'console_scripts': ['fprettify = fprettify:main']},
      py_modules=['fprettify',
          'formatting.normalizeFortranFile',
          'formatting.reformatFortranFile',
          'formatting.replacer',
          'formatting.selftest']
     )
