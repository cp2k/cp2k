# This file is part of s-dftd3.
# SPDX-Identifier: LGPL-3.0-or-later
#
# s-dftd3 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# s-dftd3 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

cc = meson.get_compiler('c')

pymod = import('python')
python = pymod.find_installation(
  get_option('python_version'),
  modules: [
    'cffi',  # Needed for generating the FFI interface
    'setuptools',  # Needed for running ffi-builder.py
  ],
)
python_dep = python.dependency(required: true)

# Python's CFFI is horrible in working with preprocessor statements,
# therefore, we have to preprocess the header before passing it to the ffibuilder
sdftd3_pp = configure_file(
  command: [
    cc,
    '-I@0@'.format(
      sdftd3_dep.get_variable(
        pkgconfig: 'includedir',
        cmake: 's-dftd3_INCLUDE_DIRS',
        internal: 'includedir',
      ).split().get(0)
    ),
    '-DSDFTD3_CFFI',
    '-E',
    '@INPUT@',
  ],
  input: sdftd3_header[0],
  output: '_libdftd3.h',
  capture: true,
)

# This is the actual out-of-line API processing of the ffibuilder
sdftd3_cffi_srcs = configure_file(
  command: [python, files('..'/'ffibuilder.py'), '@INPUT@', '@BASENAME@'],
  input: sdftd3_pp,
  output: '@BASENAME@.c',
)

# Actual generation of the Python extension, since the shared_module does not work
# well with dependency objects, we will trick it by linking a whole static lib
sdftd3_pyext = python.extension_module(
  '_libdftd3',
  link_whole: static_library(
    '_libdftd3',
    sdftd3_cffi_srcs,
    dependencies: [sdftd3_dep, python_dep],
  ),
  dependencies: [sdftd3_dep, python_dep],
  install: install,
  subdir: 'dftd3',
)

pysrcs = files(
  '__init__.py',
  'ase.py',
  'library.py',
  'interface.py',
  'parameters.py',
  'pyscf.py',
  'qcschema.py',
  'test_ase.py',
  'test_library.py',
  'test_interface.py',
  'test_parameters.py',
  'test_pyscf.py',
  'test_qcschema.py',
  'py.typed',
)
fs = import('fs')
if fs.exists('parameters.toml')
  pysrcs += files('parameters.toml')
endif

if install
  python.install_sources(
    pysrcs,
    subdir: 'dftd3',
    pure: false,
  )
endif
