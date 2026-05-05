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

if(WITH_API)
  add_subdirectory("api")
endif()
add_subdirectory("acp")
add_subdirectory("basis")
add_subdirectory("blas")
add_subdirectory("ceh")
add_subdirectory("classical")
add_subdirectory("container")
add_subdirectory("context")
add_subdirectory("coulomb")
add_subdirectory("disp")
add_subdirectory("data")
add_subdirectory("exchange")
add_subdirectory("external")
add_subdirectory("fit")
add_subdirectory("integral")
add_subdirectory("io")
add_subdirectory("lapack")
add_subdirectory("mesh")
add_subdirectory("output")
add_subdirectory("param")
add_subdirectory("post_processing")
add_subdirectory("repulsion")
add_subdirectory("scf")
add_subdirectory("solvation")
add_subdirectory("utils")
add_subdirectory("wavefunction")
add_subdirectory("xtb")

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/acp.f90"
  "${dir}/adjlist.f90"
  "${dir}/basis.f90"
  "${dir}/blas.f90"
  "${dir}/container.f90"
  "${dir}/context.f90"
  "${dir}/coulomb.f90"
  "${dir}/cutoff.f90"
  "${dir}/data.f90"
  "${dir}/double_dictionary.f90"
  "${dir}/disp.f90"
  "${dir}/exchange.f90"
  "${dir}/lapack.f90"
  "${dir}/os.F90"
  "${dir}/param.f90"
  "${dir}/repulsion.f90"
  "${dir}/results.f90"
  "${dir}/scf.f90"
  "${dir}/solvation.f90"
  "${dir}/spin.f90"
  "${dir}/timer.f90"
  "${dir}/toml.f90"
  "${dir}/utils.f90"
  "${dir}/version.f90"
  "${dir}/wavefunction.f90"
  "${dir}/wignerseitz.f90"
  "${dir}/xtb.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
