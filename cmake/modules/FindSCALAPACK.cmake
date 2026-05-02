#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(SCALAPACK "SCALAPACK")

# check if we have mkl as blas library or not and pick the scalapack from mkl
# distro if found
if(NOT CP2K_CONFIG_PACKAGE)
  if(CP2K_SCALAPACK_VENDOR MATCHES "MKL|auto"
     AND TARGET cp2k::BLAS::MKL::scalapack_link)
    # we have mkl check for the different mkl target
    get_target_property(
      CP2K_SCALAPACK_LINK_LIBRARIES cp2k::BLAS::MKL::scalapack_link
      INTERFACE_LINK_LIBRARIES)
    set(CP2K_SCALAPACK_FOUND yes)
  elseif(CP2K_SCALAPACK_VENDOR MATCHES "SCI|auto"
         AND TARGET cp2k::BLAS::SCI::scalapack_link)
    get_target_property(
      CP2K_SCALAPACK_LINK_LIBRARIES cp2k::BLAS::SCI::scalapack_link
      INTERFACE_LINK_LIBRARIES)
    set(CP2K_SCALAPACK_FOUND yes)
  elseif(CP2K_SCALAPACK_VENDOR MATCHES "NVPL")
    if(CP2K_BLAS_INTERFACE MATCHES "64bits")
      set(_nvpl_int_type "ilp64")
    else()
      set(_nvpl_int_type "lp64")
    endif()

    set(_nvpl_mpi_type "${CP2K_NVPL_SCALAPACK_MPI}")
    if(_nvpl_mpi_type STREQUAL "auto")
      if(MPI_Fortran_LIBRARY_VERSION_STRING MATCHES "Open MPI[^0-9]*([0-9]+)")
        set(_nvpl_mpi_type "openmpi${CMAKE_MATCH_1}")
      elseif(MPI_Fortran_LIBRARY_VERSION_STRING MATCHES "MPICH|HYDRA")
        set(_nvpl_mpi_type "mpich")
      else()
        message(
          FATAL_ERROR
            "Could not determine the NVPL BLACS MPI interface. Set "
            "CP2K_NVPL_SCALAPACK_MPI to mpich, openmpi3, openmpi4, or openmpi5."
        )
      endif()
    endif()

    find_package(nvpl REQUIRED COMPONENTS scalapack)
    set(_nvpl_scalapack_target "nvpl::scalapack_${_nvpl_int_type}")
    set(_nvpl_blacs_target "nvpl::blacs_${_nvpl_int_type}_${_nvpl_mpi_type}")
    if(NOT TARGET "${_nvpl_scalapack_target}")
      message(
        FATAL_ERROR "NVPL ScaLAPACK target ${_nvpl_scalapack_target} not found")
    endif()
    if(NOT TARGET "${_nvpl_blacs_target}")
      message(FATAL_ERROR "NVPL BLACS target ${_nvpl_blacs_target} not found")
    endif()
    set(CP2K_SCALAPACK_LINK_LIBRARIES
        "${_nvpl_scalapack_target};${_nvpl_blacs_target}")
    set(CP2K_SCALAPACK_FOUND yes)
  else() # if(CP2K_SCALAPACK_VENDOR MATCHES "GENERIC|auto")
    if(TARGET cp2k::BLAS::MKL::scalapack_link)
      message(
        WARNING
          "-----------------------------------------------------------------"
          "-                  FindScalapack warning                        -"
          "-----------------------------------------------------------------"
          "\n"
          "You may want to use mkl implementation of scalapack. To do this\n"
          "add -DCP2K_SCALAPACK_VENDOR=MKL to the cmake command line.\n")
    elseif(TARGET cp2k::BLAS::SCI::scalapack_link)
      message(
        WARNING
          "-----------------------------------------------------------------"
          "-                  FindScalapack warning                        -"
          "-----------------------------------------------------------------"
          "\n"
          "You may want to use Cray implementation of scalapack. To do this\n"
          "add -DCP2K_SCALAPACK_VENDOR=SCI to the cmake command line\n\n")
    endif()

    # try to detect location with pkgconfig
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(CP2K_SCALAPACK IMPORTED_TARGET GLOBAL "scalapack")
    endif()

    # this should be enough for detecting scalapack compiled by hand. If
    # scalapack is vendor specific then we sahould have a target blas::scalapack
    # available. it removes the problem of modifying too many files when we add
    # a vendor specific blas/lapack/scalapack implementation
    if(CP2K_USE_MPI)
      if(NOT CP2K_SCALAPACK_FOUND)
        if("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI")
          cp2k_find_libraries(SCALAPACK "scalapack-openmpi")
        else()
          cp2k_find_libraries(SCALAPACK "scalapack-mpich")
        endif()
      endif()
    endif()
    if(NOT CP2K_SCALAPACK_FOUND)
      cp2k_find_libraries(SCALAPACK "scalapack")
    endif()
  endif()
endif()

# cleanup list (regularly contains empty items)
list(FILTER CP2K_SCALAPACK_LINK_LIBRARIES EXCLUDE REGEX "^$")

# check if found
find_package_handle_standard_args(SCALAPACK
                                  REQUIRED_VARS CP2K_SCALAPACK_LINK_LIBRARIES)
# prevent clutter in cache

# add target to link against
if(NOT TARGET cp2k::SCALAPACK::scalapack)
  add_library(cp2k::SCALAPACK::scalapack INTERFACE IMPORTED)
endif()

set_property(TARGET cp2k::SCALAPACK::scalapack
             PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_SCALAPACK_LINK_LIBRARIES})
mark_as_advanced(CP2K_SCALAPACK_LINK_LIBRARIES)
