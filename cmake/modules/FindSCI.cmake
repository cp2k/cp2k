#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

# set paths to look for library from ROOT variables.If new policy is set,
# find_library() automatically uses them.
include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(LIBSCI "SCI")

# we might need to change the logic a little here since the cp2k_find_library
# function expect to have CP2K_package_PREFIX set.

set(CP2K_LIBSCI_MP_ROOT "${CP2K_LIBSCI_ROOT}")
set(CP2K_LIBSCI_MPI_ROOT "${CP2K_LIBSCI_ROOT}")
set(CP2K_LIBSCI_MPI_MP_ROOT "${CP2K_LIBSCI_ROOT}")

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
  if(NOT TARGET cp2k::BLAS::SCI::sci)
    add_library(cp2k::BLAS::SCI::sci INTERFACE IMPORTED)
    add_library(cp2k::BLAS::SCI::sci_mpi INTERFACE IMPORTED)
    add_library(cp2k::BLAS::SCI::sci_mp INTERFACE IMPORTED)
    add_library(cp2k::BLAS::SCI::sci_mpi_mp INTERFACE IMPORTED)
    add_library(cp2k::BLAS::SCI::scalapack_link INTERFACE IMPORTED)
    add_library(cp2k::BLAS::SCI::blas INTERFACE IMPORTED)

    if(CP2K_LIBSCI_INCLUDE_DIRS)
      set_property(
        TARGET cp2k::BLAS::SCI::sci PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                             "${CP2K_LIBSCI_INCLUDE_DIRS}")
      set_property(
        TARGET cp2k::BLAS::SCI::sci_mp PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                "${CP2K_LIBSCI_INCLUDE_DIRS}")
      set_property(
        TARGET cp2k::BLAS::SCI::sci_mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                 "${CP2K_LIBSCI_INCLUDE_DIRS}")
      set_property(
        TARGET cp2k::BLAS::SCI::sci_mpi_mp
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBSCI_INCLUDE_DIRS}")
    endif()

    set_property(
      TARGET cp2k::BLAS::SCI::sci PROPERTY INTERFACE_LINK_LIBRARIES
                                           ${CP2K_LIBSCI_LINK_LIBRARIES})
    set_property(
      TARGET cp2k::BLAS::SCI::sci_mp PROPERTY INTERFACE_LINK_LIBRARIES
                                              ${CP2K_LIBSCI_MP_LINK_LIBRARIES})
    set_property(
      TARGET cp2k::BLAS::SCI::sci_mpi
      PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_LIBSCI_MPI_LINK_LIBRARIES}
               cp2k::BLAS::SCI::sci)
    set_property(
      TARGET cp2k::BLAS::SCI::sci_mpi_mp
      PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_LIBSCI_MPI_MP_LINK_LIBRARIES}
               cp2k::BLAS::SCI::sci_mp)
    set_property(
      TARGET cp2k::BLAS::SCI::scalapack_link
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBSCI_INCLUDE_DIRS}")
  endif()

  if(CP2K_BLAS_THREADING MATCHES "sequential")
    set_property(TARGET cp2k::BLAS::SCI::blas PROPERTY INTERFACE_LINK_LIBRARIES
                                                       cp2k::BLAS::SCI::sci)
    set_property(TARGET cp2k::BLAS::SCI::scalapack_link
                 PROPERTY INTERFACE_LINK_LIBRARIES cp2k::BLAS::SCI::sci_mpi)
  else()
    set_property(TARGET cp2k::BLAS::SCI::blas PROPERTY INTERFACE_LINK_LIBRARIES
                                                       cp2k::BLAS::SCI::sci_mp)
    set_property(TARGET cp2k::BLAS::SCI::scalapack_link
                 PROPERTY INTERFACE_LINK_LIBRARIES cp2k::BLAS::SCI::sci_mpi_mp)
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
