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

add_subdirectory("damping")
add_subdirectory("data")
add_subdirectory("integrator")
add_subdirectory("model")
add_subdirectory("param")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/blas.F90"
  "${dir}/cache.f90"
  "${dir}/cutoff.f90"
  "${dir}/damping.f90"
  "${dir}/data.f90"
  "${dir}/disp.f90"
  "${dir}/integrator.f90"
  "${dir}/model.f90"
  "${dir}/ncoord.f90"
  "${dir}/numdiff.f90"
  "${dir}/output.f90"
  "${dir}/param.f90"
  "${dir}/utils.f90"
  "${dir}/version.f90"
)
if(WITH_API)
  list(APPEND srcs "${dir}/api.f90")
endif()

set(srcs "${srcs}" PARENT_SCOPE)
