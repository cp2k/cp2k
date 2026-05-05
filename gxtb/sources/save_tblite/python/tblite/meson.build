# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

cc = meson.get_compiler('c')

pymod = import('python')
python = pymod.find_installation(
  get_option('python_version'),
  modules: [
    'cffi',
    'setuptools',
  ],
)
python_dep = python.dependency(required: true)

# Python's CFFI is horrible in working with preprocessor statements,
# therefore, we have to preprocess the header before passing it to the ffibuilder
tblite_pp = configure_file(
  command: [
    cc,
    '-I@0@'.format(
      tblite_dep.get_variable(
        pkgconfig: 'includedir',
        cmake: 'tblite_INCLUDE_DIRS',
        internal: 'includedir',
      ).split().get(0)
    ),
    '-DTBLITE_CFFI',
    '-E',
    '@INPUT@',
  ],
  input: tblite_header,
  output: '_libtblite.h',
  capture: true,
)

# This is the actual out-of-line API processing of the ffibuilder
tblite_cffi_srcs = configure_file(
  command: [python, files('..'/'ffi-builder.py'), '@INPUT@', '@BASENAME@'],
  input: tblite_pp,
  output: '@BASENAME@.c',
)

# Actual generation of the Python extension, since the shared_module does not work
# well with dependency objects, we will trick it by linking a whole static lib
tblite_pyext = python.extension_module(
  '_libtblite',
  link_whole: static_library(
    '_libtblite',
    tblite_cffi_srcs,
    dependencies: [tblite_dep, python_dep],
  ),
  dependencies: [tblite_dep, python_dep],
  install: install,
  subdir: 'tblite',
  link_language: 'c',
)

if install
  python.install_sources(
    files(
      '__init__.py',
      'ase.py',
      'exceptions.py',
      'interface.py',
      'library.py',
      'qcschema.py',
      'test_ase.py',
      'test_interface.py',
      'test_library.py',
      'test_qcschema.py',
      'py.typed',
    ),
    subdir: 'tblite',
    pure: false,
  )
endif

