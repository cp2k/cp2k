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

# Standalone build for Python-API of DFT-D4, requires a dftd4 library
# in the PKG_CONFIG_PATH environment variable to work.
project(
  'dftd3',
  'c',
  version: '1.3.2',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.55,!=1.8.0',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
install = true

sdftd3_dep = dependency(
  's-dftd3',
  version: '>=@0@'.format(meson.project_version()),
  fallback: ['s-dftd3', 'sdftd3_dep'],
  default_options: [
    'default_library=static',
    'api=true',
    'python=false',
  ],
)
sdftd3_header = files('include'/'_dftd3.h')

subdir('dftd3')
