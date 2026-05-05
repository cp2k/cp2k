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

sdftd3_exe = executable(
  meson.project_name(),
  sources: files(
    'argument.f90',
    'cli.f90',
    'driver.f90',
    'help.f90',
    'main.f90',
    'toml.f90',
  ),
  dependencies: [sdftd3_dep, exe_deps],
  install: install,
)

if os != 'windows'
  test('version', sdftd3_exe, args: '--version', suite: 'app')
  test('help', sdftd3_exe, args: '--help', suite: 'app')
  test('noargs', sdftd3_exe, should_fail: true, suite: 'app')
  test('example', sdftd3_exe, args: files('example.xyz'), suite: 'app')
  test('run-example', sdftd3_exe, args: ['run', '--', files('example.xyz')], suite: 'app')
  test('param-all', sdftd3_exe, args: ['param', sdftd3_parameters], suite: 'app')
  test('param-select', sdftd3_exe, args: ['param', sdftd3_parameters, 'b2plyp'], suite: 'app')
endif
