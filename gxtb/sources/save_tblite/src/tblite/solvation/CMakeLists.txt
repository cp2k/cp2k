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

add_subdirectory("data")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/alpb.f90"
  "${dir}/born.f90"
  "${dir}/cds.f90"
  "${dir}/cm5.f90"
  "${dir}/cpcm.f90"
  "${dir}/cpcm_dd.f90"
  "${dir}/data.f90"
  "${dir}/input.f90"
  "${dir}/surface.f90"
  "${dir}/shift.f90"
  "${dir}/type.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
