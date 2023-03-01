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

# check for blas first. Most of the vendor libraries bundle lapack and blas in
# the same library. If so the FindBlas.cmake module will contain this
# information already and the information will be included in the blas target
#
# This solution might not good enough though.

find_package(PkgConfig)
find_package(Blas REQUIRED)

if(CP2K_BLAS_FOUND)
  # LAPACK in the Intel MKL 10+ library?
  if(CP2K_BLAS_VENDOR MATCHES "MKL|OpenBLAS|Armpl|SCI|FlexiBLAS|NVHPC")
    # we just need to create the interface that's all
    set(CP2K_LAPACK_FOUND TRUE)
    get_target_property(CP2K_LAPACK_INCLUDE_DIRS CP2K::BLAS::blas
                        INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(CP2K_LAPACK_LIBRARIES CP2K::BLAS::blas
                        INTERFACE_LINK_LIBRARIES)
  else()
    # we might get lucky to find a pkgconfig package for lapack (fedora provides
    # one for instance)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(CP2K_LAPACK lapack)
    endif()

    if(NOT CP2K_LAPACK_FOUND)
      find_library(
        CP2K_LAPACK_LIBRARIES
        NAMES "lapack" "lapack64"
        PATH_SUFFIXES "openblas" "openblas64" "openblas-pthread"
                      "openblas-openmp" "lib" "lib64"
        NO_DEFAULT_PATH)
    endif()
  endif()
endif()

# check if found
find_package_handle_standard_args(Lapack REQUIRED_VARS CP2K_LAPACK_LIBRARIES)

if(NOT TARGET CP2K::LAPACK::lapack)
  add_library(CP2K::LAPACK::lapack INTERFACE IMPORTED)
  add_library(CP2K::LAPACK::LAPACK ALIAS CP2K::LAPACK::lapack)
endif()
set_property(TARGET CP2K::LAPACK::lapack PROPERTY INTERFACE_LINK_LIBRARIES
                                                  ${CP2K_LAPACK_LIBRARIES})
if(CP2K_LAPACK_INCLUDE_DIRS)
  set_property(
    TARGET CP2K::LAPACK::lapack PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                         ${CP2K_LAPACK_INCLUDE_DIRS})
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_LAPACK_FOUND CP2K_LAPACK_LIBRARIES
                 CP2K_LAPACK_INCLUDE_DIRS)
