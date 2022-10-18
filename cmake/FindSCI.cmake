# Copyright (c) 2019- ETH Zurich
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# set paths to look for library from ROOT variables.If new policy is set,
# find_library() automatically uses them.
include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(LIBSCI)

# we might need to change the logic a little here since the cp2k_find_library
# function expect to have CP2K_package_PREFIX set.

set(CP2K_LIBSCI_MP_PREFIX "${CP2K_LIBSCI_PREFIX}")
set(CP2K_LIBSCI_MPI_PREFIX "${CP2K_LIBSCI_PREFIX}")
set(CP2K_LIBSCI_MPI_MP_PREFIX "${CP2K_LIBSCI_PREFIX}")

set(_sci_lib "sci_gnu")

if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set(_sci_lib "sci_intel")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(_sci_lib "sci_cray")
endif()

cp2k_find_libraries("LIBSCI" "${_sci_lib}")
cp2k_find_libraries("LIBSCI_MP" "${_sci_lib}_mp")
cp2k_find_libraries("LIBSCI_MPI" "${_sci_lib}_mpi")
cp2k_find_libraries("LIBSCI_MPI_MP" "${_sci_lib}_mpi_mp")
cp2k_include_dirs(LIBSCI "cblas.h")

# check if found
find_package_handle_standard_args(SCI REQUIRED_VARS CP2K_LIBSCI_INCLUDE_DIRS
                                                    CP2K_LIBSCI_LINK_LIBRARIES)

# add target to link against
if(CP2K_LIBSCI_FOUND)
  if(NOT TARGET CP2K_SCI::sci)
    add_library(CP2K_SCI::sci INTERFACE IMPORTED)
    add_library(CP2K_SCI::sci_mpi INTERFACE IMPORTED)
    add_library(CP2K_SCI::sci_mp INTERFACE IMPORTED)
    add_library(CP2K_SCI::sci_mpi_mp INTERFACE IMPORTED)
    add_library(CP2K_SCI::scalapack_link INTERFACE IMPORTED)
    add_library(CP2K_SCI::blas INTERFACE IMPORTED)

    if(CP2K_LIBSCI_INCLUDE_DIRS)
      set_property(TARGET CP2K_SCI::sci PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                 "${CP2K_LIBSCI_INCLUDE_DIRS}")
      set_property(
        TARGET CP2K_SCI::sci_mp PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                         "${CP2K_LIBSCI_INCLUDE_DIRS}")
      set_property(
        TARGET CP2K_SCI::sci_mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                          "${CP2K_LIBSCI_INCLUDE_DIRS}")
      set_property(
        TARGET CP2K_SCI::sci_mpi_mp PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                             "${CP2K_LIBSCI_INCLUDE_DIRS}")
    endif()

    set_property(TARGET CP2K_SCI::sci PROPERTY INTERFACE_LINK_LIBRARIES
                                               ${CP2K_LIBSCI_LINK_LIBRARIES})
    set_property(
      TARGET CP2K_SCI::sci_mp PROPERTY INTERFACE_LINK_LIBRARIES
                                       ${CP2K_LIBSCI_MP_LINK_LIBRARIES})
    set_property(
      TARGET CP2K_SCI::sci_mpi
      PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_LIBSCI_MPI_LINK_LIBRARIES}
               CP2K_SCI::sci)
    set_property(
      TARGET CP2K_SCI::sci_mpi_mp
      PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_LIBSCI_MPI_MP_LINK_LIBRARIES}
               CP2K_SCI::sci_mp)
    set_property(
      TARGET CP2K_SCI::scalapack_link PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                               "${CP2K_LIBSCI_INCLUDE_DIRS}")
  endif()

  if(CP2K_BLAS_THREADING MATCHES "sequential")
    set_property(TARGET CP2K_SCI::blas PROPERTY INTERFACE_LINK_LIBRARIES
                                                CP2K_SCI::sci)
    set_property(TARGET CP2K_SCI::scalapack_link
                 PROPERTY INTERFACE_LINK_LIBRARIES CP2K_SCI::sci_mpi)
  else()
    set_property(TARGET CP2K_SCI::blas PROPERTY INTERFACE_LINK_LIBRARIES
                                                CP2K_SCI::sci_mp)
    set_property(TARGET CP2K_SCI::scalapack_link
                 PROPERTY INTERFACE_LINK_LIBRARIES CP2K_SCI::sci_mpi_mp)
  endif()

  set(CP2K_BLAS_VENDOR "SCI")

endif()

# prevent clutter in cache
mark_as_advanced(
  CP2K_LIBSCI_FOUND
  CP2K_BLAS_VENDOR
  CP2K_LIBSCI_LINK_LIBRARIES
  CP2K_LIBSCI_MP_LINK_LIBRARIES
  CP2K_LIBSCI_MPI_LINK_LIBRARIES
  CP2K_LIBSCI_MPI_MP_LINK_LIBRARIES
  CP2K_LIBSCI_INCLUDE_DIRS)
