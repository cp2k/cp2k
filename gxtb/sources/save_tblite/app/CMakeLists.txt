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

add_executable(
  "${PROJECT_NAME}-exe"
  "argument.f90"
  "cli.f90"
  "cli_help.f90"
  "driver.f90"
  "driver_fit.f90"
  "driver_param.f90"
  "driver_run.f90"
  "driver_guess.f90"
  "driver_tagdiff.f90"
  "features.F90"
  "main.f90"
)
set_target_properties(
  "${PROJECT_NAME}-exe"
  PROPERTIES
  OUTPUT_NAME "${PROJECT_NAME}"
)
target_compile_options(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "$<$<Fortran_COMPILER_ID:GNU>:-fall-intrinsics>"
)
target_link_libraries(
  "${PROJECT_NAME}-exe"
  PRIVATE
  "${PROJECT_NAME}-lib"
)

install(
  TARGETS
  "${PROJECT_NAME}-exe"
  DESTINATION
  "${CMAKE_INSTALL_BINDIR}"
)
