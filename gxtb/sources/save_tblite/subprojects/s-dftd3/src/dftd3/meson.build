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

subdir('damping')
subdir('data')
subdir('gcp')

srcs += files(
  'citation.f90',
  'cutoff.f90',
  'damping.f90',
  'data.f90',
  'disp.f90',
  'gcp.f90',
  'model.f90',
  'ncoord.f90',
  'output.f90',
  'param.f90',
  'reference.f90',
  'utils.f90',
  'version.f90',
)
if get_option('api') or get_option('python')
  srcs += files('api.f90')
endif
