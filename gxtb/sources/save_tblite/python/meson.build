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

project(
  'tblite',
  'c',
  version: '0.5.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.57.2,!=1.8.0',
  default_options: [
    'buildtype=debugoptimized',
  ],
)
install = true

tblite_dep = dependency(
  'tblite',
  version: '>=@0@'.format(meson.project_version()),
  fallback: ['tblite', 'tblite_dep'],
  default_options: [
    'default_library=static',
    'api=true',
    'python=false',
  ],
)
tblite_header = files('include'/'_tblite.h')

subdir('tblite')
