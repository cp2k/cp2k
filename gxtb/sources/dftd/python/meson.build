# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

# Standalone build for Python-API of DFT-D4, requires a dftd4 library
# in the PKG_CONFIG_PATH environment variable to work.
project(
  'dftd4',
  'c',
  version: '4.0.2',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.55,!=1.8.0',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
install = true

dftd4_dep = dependency(
  'dftd4',
  version: '>=@0@'.format(meson.project_version()),
  fallback: ['dftd4', 'dftd4_dep'],
  default_options: [
    'default_library=static',
    'api=true',
    'python=false',
  ],
)
dftd4_header = files('include'/'_dftd4.h')

subdir('dftd4')
