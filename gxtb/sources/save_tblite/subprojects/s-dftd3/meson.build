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

project(
  's-dftd3',
  'fortran',
  version: '1.3.2',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.55,!=1.8.0',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
    'c_std=c11',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: get_option('api') or get_option('python'))

# General configuration information
inc_dirs = []
lib_deps = []
exe_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
sdftd3_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
sdftd3_inc = [sdftd3_lib.private_dir_include(), include_directories('include')]
sdftd3_dep = declare_dependency(
  link_with: sdftd3_lib,
  include_directories: sdftd3_inc,
  dependencies: lib_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)

sdftd3_parameters = files(
  'assets/parameters.toml',
)

# Add applications
subdir('app')

# Package the license files
sdftd3_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

sdftd3_header = files(
  'include/dftd3.h',
  'include/s-dftd3.h',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    sdftd3_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name(),
  )

  if get_option('api')
    install_headers(
      sdftd3_header,
    )
  endif

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    sdftd3_lib,
    description: 'Simple reimplementation of the DFT-D3 dispersion model',
    subdirs: ['', module_id],
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man/s-dftd3.1.adoc'),
        output: '@BASENAME@',
      )
    )
  endif

  install_data(
    sdftd3_parameters,
    install_dir: get_option('datadir')/meson.project_name(),
  )
endif

# Add the testsuite
subdir('test')

# Build Python extension module
if get_option('python')
  subdir('python'/'dftd3')
endif
