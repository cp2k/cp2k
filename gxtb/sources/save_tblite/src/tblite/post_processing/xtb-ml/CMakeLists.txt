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
  "${dir}/atomic_frontier_orbitals.f90"
  "${dir}/convolution.f90"
  "${dir}/density.f90"
  "${dir}/energy.f90"
  "${dir}/features.f90"
  "${dir}/geometry.f90"
  "${dir}/orbital_energy.f90"
  "${dir}/xtbml.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
