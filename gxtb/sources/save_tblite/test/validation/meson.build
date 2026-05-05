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

validation_tester = find_program(files('tester.py'))

run_element = [
  '01_h3+',
  '02_heada',
  '03_etlicl-',
  '04_nebeo',
  '05_b3o3me3',
  '06_bic221',
  '07_bipyrrol',
  '08_ketene',
  '09_fccf',
  '10_nenacl',
  '11_crown6+',
  '12_mgcl2h2o2',
  '13_al3n3h6',
  '14_sif5-',
  '15_ylid',
  '16_fshsf2',
  '17_c6cl6',
  '18_nebz',
  '19_knapht+',
  '20_cacp2',
]

foreach m : ['gfn2', 'gfn1', 'ipea1']
  foreach t : run_element
    benchmark(
      '@0@-@1@'.format(m, t),
      validation_tester,
      args: [
        tblite_exe,
        files(t/'@0@-xtb.json'.format(m)),
        files(t/'coord'),
        '--no-restart',
        '--method', m,
        '--grad', '@0@-@1@.txt'.format(m, t),
        t.contains('+') ? ['--charge', '+1'] : [],
        t.contains('-') ? ['--charge', '-1'] : [],
      ],
      suite: 'calc',
    )
  endforeach
endforeach
