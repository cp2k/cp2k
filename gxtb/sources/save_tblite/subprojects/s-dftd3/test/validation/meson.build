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

app_tester = find_program(meson.current_source_dir() / 'tester.py')
tests = [
  '01-energy-d3-bj',
  '02-energy-d3-zero',
  '03-energy-d3-bjm',
  '04-energy-d3-zerom',
  '05-energy-d3-op',
  '06-gradient-d3-bj',
  '07-gradient-d3-zero',
  '08-gradient-d3-bjm',
  '09-gradient-d3-zerom',
  '10-gradient-d3-op',
  '11-pairwise-d3-bj',
  '12-pairwise-d3-zero',
  '13-pairwise-d3-bjm',
  '14-pairwise-d3-zerom',
  '15-pairwise-d3-op',
  '16-db-d3-bj',
  '17-db-d3-zero',
  '18-db-d3-bjm',
  '19-db-d3-zerom',
  '20-db-d3-op',
  '21-gcp-pbeh3c',
  '22-gcp-b973c',
  '23-gcp-hf3c',
  '24-gcp-hf-deftzvp',
  '25-gcp-b3lyp-631gd',
]

foreach t : tests
  benchmark(
    t,
    app_tester,
    args: [
      sdftd3_exe,
      files('@0@.json'.format(t)),
      files('@0@.resp'.format(t)),
    ],
    suite: 'app',
  )
endforeach
