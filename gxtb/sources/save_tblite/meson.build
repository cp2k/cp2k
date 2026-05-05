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
  'fortran',
  version: '0.5.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.57.2,!=1.8.0',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: get_option('api') or get_option('python'), native: false)

# General configuration information
lib_deps = []
inc_dirs = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
tblite_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
  link_language: 'fortran',
)

# Export dependency for other projects and test suite
tblite_inc = [include_directories('include'), tblite_lib.private_dir_include()]
tblite_dep = declare_dependency(
  link_with: tblite_lib,
  include_directories: tblite_inc,
  dependencies: lib_deps,
  variables: {'includedir': meson.current_source_dir() / 'include'},
)

# Add executable targets
subdir('app')

# Package the license files
tblite_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

tblite_header = files(
  'include/tblite.h',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    tblite_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  install_headers(
    tblite_header,
  )
  install_subdir(
    'include/tblite',
    install_dir: get_option('includedir'),
  )

  module_id = meson.project_name() / fc_id + '-' + fc.version()
  meson.add_install_script(
    find_program(files('config'/'install-mod.py')),
    get_option('includedir') / module_id,
  )

  pkg = import('pkgconfig')
  pkg.generate(
    tblite_lib,
    description: 'Light-weight tight-binding framework',
    subdirs: ['', module_id],
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    man_srcs = files(
      'man/tblite.1.adoc',
      'man/tblite-fit.1.adoc',
      'man/tblite-run.1.adoc',
      'man/tblite-guess.1.adoc',
      'man/tblite-param.1.adoc',
      'man/tblite-tagdiff.1.adoc',
      'man/tblite-tag.5.adoc',
      'man/tblite-solvents.7.adoc',
    )
    foreach man : man_srcs
      install_man(
        configure_file(
          command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
          input: man,
          output: '@BASENAME@',
        )
      )
    endforeach
  endif
endif

# add the testsuite
subdir('test')

if get_option('python')
  subdir('python/tblite')
endif
