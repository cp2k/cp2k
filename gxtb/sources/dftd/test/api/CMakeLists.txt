# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

add_executable(
  "${PROJECT_NAME}-api-tester"
  "example.c"
)
target_include_directories(
  "${PROJECT_NAME}-api-tester"
  PRIVATE
  "${PROJECT_SOURCE_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-api-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
)
add_test(NAME "api" COMMAND "${PROJECT_NAME}-api-tester")
