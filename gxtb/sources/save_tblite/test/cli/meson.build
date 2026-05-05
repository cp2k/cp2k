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

# Test help pages for each subprogram
test('help', tblite_exe, args: '--help', suite: 'cli')
foreach subprog: ['fit', 'param', 'run', 'guess', 'tagdiff']
  test('@0@-help'.format(subprog), tblite_exe, args: [subprog, '--help'], suite: 'cli')
endforeach

# Test version for each subprogram
test('version', tblite_exe, args: '--version', suite: 'cli')
foreach subprog: ['fit', 'param', 'run', 'guess', 'tagdiff']
  test('@0@-version'.format(subprog), tblite_exe, args: [subprog, '--version'], suite: 'cli')
endforeach

# Test that the program fails when no arguments are given
test('noargs', tblite_exe, should_fail: true, suite: 'cli')
foreach subprog: ['fit', 'param', 'run', 'guess', 'tagdiff']
  test('@0@-noargs'.format(subprog), tblite_exe, args: subprog, should_fail: true, suite: 'cli')
endforeach

app_tester = find_program(files('tester.py'))

run_testcases = [
  '01a-gfn1xtb',
  '01b-gfn2xtb',
  '01c-ipea1xtb',
  '01d-gxtb',
  '01e-gfn1xtb-solv',
  '01f-gfn2xtb-solv',
  '02a-gfn1xtb',
  '02b-gfn2xtb',
  '02c-ipea1xtb',
  '02d-gxtb',
  '02e-gfn1xtb-solv',
  '02f-gfn2xtb-solv',
  '03a-gfn1xtb',
  '03b-gfn2xtb',
  '03c-ipea1xtb',
  '03d-gxtb',
  '03e-gfn1xtb-solv',
  '03f-gfn2xtb-solv',
  '04a-gfn1xtb',
  '04b-gfn2xtb',
  '04c-ipea1xtb',
  '04d-gxtb',
  '04e-gfn1xtb-solv',
  '04f-gfn2xtb-solv',
  '05a-gfn1xtb',
  '05b-gfn2xtb',
  '05c-ipea1xtb',
  '05d-gxtb',
  '05e-gfn1xtb-solv',
  '05f-gfn2xtb-solv',
  '06a-gfn1xtb',
  '07f-ceh',
  '08d-gxtb',
  '09b-gfn2xtb',
  '09d-gxtb',
]
foreach case: run_testcases
  test(
    case,
    app_tester,
    args: [
      tblite_exe,
      files('@0@.json'.format(case)),
      files('@0@.resp'.format(case)),
    ],
    suite: 'cli',
  )
endforeach

test(
  'tagdiff-example',
  tblite_exe,
  args: [
    'tagdiff',
    files('test'/'results.tag'),
    files('test'/'reference.tag'),
  ],
  suite: 'cli',
)

test(
  'fit-example',
  tblite_exe,
  args: [
    'fit', '--dry-run',
    files('test'/'simple-xtb.toml'),
    files('test'/'fit-input.toml'),
  ],
  suite: 'cli',
)

foreach m : ['gfn2', 'gfn1', 'ipea1']
  test(
    'param-@0@'.format(m),
    tblite_exe,
    args: [
      'param',
      '--method', m,
      '--output', '@0@-xtb.toml'.format(m),
    ],
    suite: 'cli',
  )
endforeach

test(
  'post-processing-str', 
  tblite_exe, 
  args: [
    'run', 
    '--no-restart',
    '--post-processing', 
    'molmom', 
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)

test(
  'post-processing-toml', 
  tblite_exe, 
  args: [
    'run',
    '--no-restart',
    '--post-processing', 
    files('test'/'molmom.toml'), 
    '-v',
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)

test(
  'post-processing-param', 
  tblite_exe, 
  args: [
    'run', 
    '--no-restart',
    '--param', 
    files('test'/'gfn2-molmom.toml'), 
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)

test(
  'post-processing-other-file', 
  tblite_exe, 
  should_fail:true, 
  args: [
    'run', 
    '--no-restart',
    '--post-processing', 
    files('test'/'reference.tag'), 
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)

test(
  'xtbml-cli-arg', 
  tblite_exe, 
  should_fail:false, 
  args: [
    'run', 
    '--no-restart',
    '--post-processing', 
    'xtbml', 
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)

test(
  'xtbml-xyz-cli-arg', 
  tblite_exe, 
  should_fail:false, 
  args: [
    'run', 
    '--no-restart',
    '--post-processing', 
    'xtbml-xyz', 
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)

test(
  'xtbml-xyz-cli-toml', 
  tblite_exe, 
  should_fail:false, 
  args: [
    'run', 
    '--no-restart',
    '--post-processing', 
    files('test'/'xtbml.toml'), 
    files('07-mb33.xyz')
  ], 
  suite: 'cli',
)