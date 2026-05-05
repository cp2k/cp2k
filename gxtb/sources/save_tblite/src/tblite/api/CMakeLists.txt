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

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/calculator.f90"
  "${dir}/container.f90"
  "${dir}/context.f90"
  "${dir}/double_dictionary.f90"
  "${dir}/error.f90"
  "${dir}/param.f90"
  "${dir}/result.f90"
  "${dir}/solvation.f90"
  "${dir}/structure.f90"
  "${dir}/table.f90"
  "${dir}/utils.f90"
  "${dir}/version.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
