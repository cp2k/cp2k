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

dftd4_exe = executable(
  meson.project_name(),
  sources: files(
    'main.f90',
    'argument.f90',
    'cli.f90',
    'driver.f90',
    'help.f90',
  ),
  dependencies: dftd4_dep,
  install: install,
  link_language: 'fortran',
)

test('app-version', dftd4_exe, args: '--version')
test('app-help', dftd4_exe, args: '--help')
test('app-help-run', dftd4_exe, args: ['run', '--help'])
test('app-help-param', dftd4_exe, args: ['param', '--help'])
test('app-param-list', dftd4_exe, args: ['param', '--list'])
test('app-license', dftd4_exe, args: '--license')
test('app-citation', dftd4_exe, args: '--citation')
test('app-noargs', dftd4_exe, should_fail: true)

app_tester = find_program(files('tester.py'))

test(
  'app-energy',
  app_tester,
  args: [
    dftd4_exe,
    files('01-energy.json'),
    files('01-ammonia.tmol'),
    '--noedisp',
    '--func', 'tpssh'
  ],
)
test(
  'app-energy-d4s',
  app_tester,
  args: [
    dftd4_exe,
    files('01-energy-d4s.json'),
    files('01-ammonia.tmol'),
    '--noedisp',
    '--func', 'tpssh',
    '--model', 'D4S'
  ],
)
test(
  'app-gradient',
  app_tester,
  args: [
    dftd4_exe,
    files('02-gradient.json'),
    '--noedisp',
    '--func', 'scan',
    '--grad',
    '--',
    files('02-nitralin.mol'),
  ],
)
test(
  'app-gradient-d4s',
  app_tester,
  args: [
    dftd4_exe,
    files('02-gradient-d4s.json'),
    '--noedisp',
    '--func', 'scan',
    '--model', 'D4S',
    '--grad',
    '--',
    files('02-nitralin.mol'),
  ],
)
test(
  'app-properties',
  app_tester,
  args: [
    dftd4_exe,
    files('03-properties.json'),
    files('03-lenalidomid.gen'),
  ],
)
test(
  'app-pair-analysis',
  app_tester,
  args: [
    dftd4_exe,
    files('04-pair-analysis.json'),
    files('04-caffeine.xyz'),
    '--noedisp',
    '--pair-resolved',
    '--func', 'b3lyp',
  ],
)
test(
  'app-energy-d4-eeqbc',
  app_tester,
  args: [
    dftd4_exe,
    files('04-energy-d4-eeqbc.json'),
    files('04-caffeine.xyz'),
    '--noedisp',
    '--func', 'tpssh', 
    '--qmodel', 'eeqbc'
  ],
)

# r2SCAN-3c modifies s9, ga, gc, hence, requiring additional tests:
# - selection via functional name
# - check if s9, ga, gc can still be modified when name is given
# - full manual specification via command line
test(
  'app-r2scan3c',
  app_tester,
  args: [
    dftd4_exe,
    files('05-r2scan3c.json'),
    files('04-caffeine.xyz'),
    '--noedisp',
    '--func', 'r2scan-3c',
  ],
)
test(
  'app-r2scan3c-mod',
  app_tester,
  args: [
    dftd4_exe,
    files('06-r2scan3c-mod.json'),
    files('04-caffeine.xyz'),
    '--noedisp',
    '--func', 'r2scan-3c',
    '--zeta', '4.0', '5.0',
    '--s9', '1.0',
  ],
)
test(
  'app-r2scan3c-cli',
  app_tester,
  args: [
    dftd4_exe,
    files('05-r2scan3c.json'),
    files('04-caffeine.xyz'),
    '--noedisp',
    '--s6', '1.0', 
    '--s8', '0.0',
    '--a1', '0.42',
    '--a2', '5.65',
    '--s9', '2.0',
    '--rs9', '1.0',
    '--alp', '16.0',
    '--zeta', '2.0', '1.0',
  ],
)
