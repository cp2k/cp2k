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

if get_option('api') or get_option('python')
  subdir('api')
endif
subdir('acp')
subdir('basis')
subdir('blas')
subdir('classical')
subdir('container')
subdir('context')
subdir('coulomb')
subdir('ceh')
subdir('data')
subdir('disp')
subdir('exchange')
subdir('external')
subdir('fit')
subdir('integral')
subdir('io')
subdir('lapack')
subdir('mesh')
subdir('output')
subdir('param')
subdir('post_processing')
subdir('repulsion')
subdir('scf')
subdir('solvation')
subdir('utils')
subdir('wavefunction')
subdir('xtb')

srcs += files(
  'acp.f90',
  'adjlist.f90',
  'basis.f90',
  'blas.f90',
  'container.f90',
  'context.f90',
  'coulomb.f90',
  'cutoff.f90',
  'data.f90',
  'double_dictionary.f90',
  'disp.f90',
  'exchange.f90',
  'lapack.f90',
  'os.F90',
  'param.f90',
  'repulsion.f90',
  'results.f90',
  'scf.f90',
  'solvation.f90',
  'spin.f90',
  'timer.f90',
  'toml.f90',
  'utils.f90',
  'version.f90',
  'wavefunction.f90',
  'wignerseitz.f90',
  'xtb.f90',
)
