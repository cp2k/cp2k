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
  "${dir}/coulomb.f90"
  "${dir}/gfn1.f90"
  "${dir}/gfn2.f90"
  "${dir}/gxtb.f90"
  "${dir}/h0.f90"
  "${dir}/ipea1.f90"
  "${dir}/singlepoint.f90"
  "${dir}/spec.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
