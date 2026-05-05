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

exe_args = []
if fc_id == 'gcc'
  exe_args += '-fall-intrinsics'
endif

tblite_exe = executable(
  meson.project_name(),
  sources: files(
    'argument.f90',
    'cli.f90',
    'cli_help.f90',
    'driver.f90',
    'driver_fit.f90',
    'driver_param.f90',
    'driver_run.f90',
    'driver_guess.f90',
    'driver_tagdiff.f90',
    'features.F90',
    'main.f90',
  ),
  dependencies: tblite_dep,
  install: install,
  fortran_args: exe_args,
  link_language: 'fortran',
)