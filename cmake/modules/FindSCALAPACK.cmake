#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
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
  if(CP2K_SCALAPACK_VENDOR STREQUAL "GENERIC")
    if(TARGET cp2k::BLAS::MKL::scalapack_link)
      message(
        WARNING
          "-----------------------------------------------------------------"
          "-                  FindScalapack warning                        -"
          "-----------------------------------------------------------------"
          "\n"
          "You may want to use mkl implementation of scalapack. To do this\n"
          "add -DCP2K_SCALAPACK_VENDOR=MKL to the cmake command line.      \n")
    endif()

    if(TARGET cp2k::BLAS::SCI::scalapack_link)
      message(
        WARNING
          "-----------------------------------------------------------------\n"
          "-                  FindScalapack warning                        -\n"
          "-----------------------------------------------------------------\n"
          "\n"
          "You may want to use Cray implementation of scalapack. To do this \n"
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

    if(NOT CP2K_SCALAPACK_FOUND)
      cp2k_find_libraries(SCALAPACK "scalapack")
    endif()
  elseif(TARGET cp2k::BLAS::MKL::scalapack_link)
    # we have mkl check for the different mkl target
    get_target_property(
      CP2K_SCALAPACK_LINK_LIBRARIES cp2k::BLAS::MKL::scalapack_link
      INTERFACE_LINK_LIBRARIES)
    set(CP2K_SCALAPACK_FOUND yes)
  elseif(TARGET cp2k::BLAS::SCI::scalapack_link)
    # we have mkl check for the different mkl target
    get_target_property(
      CP2K_SCALAPACK_LINK_LIBRARIES cp2k::BLAS::SCI::scalapack_link
      INTERFACE_LINK_LIBRARIES)
  endif()
endif()
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
