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

add_subdirectory("damping")
add_subdirectory("data")
add_subdirectory("gcp")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/citation.f90"
  "${dir}/cutoff.f90"
  "${dir}/damping.f90"
  "${dir}/data.f90"
  "${dir}/disp.f90"
  "${dir}/gcp.f90"
  "${dir}/model.f90"
  "${dir}/ncoord.f90"
  "${dir}/output.f90"
  "${dir}/param.f90"
  "${dir}/reference.f90"
  "${dir}/utils.f90"
  "${dir}/version.f90"
)
if(WITH_API)
  list(APPEND srcs "${dir}/api.f90")
endif()

set(srcs "${srcs}" PARENT_SCOPE)
