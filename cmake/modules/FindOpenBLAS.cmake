#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(cp2k_utils)
include(FindPackageHandleStandardArgs)

find_package(PkgConfig)

cp2k_set_default_paths(OPENBLAS "OpenBLAS")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_OPENBLAS IMPORTED_TARGET GLOBAL openblas)
endif()

# try the openblas module of openblas library Maybe we are lucky it is installed
# find_package(OPENBLAS QUIET)

if(NOT CP2K_OPENBLAS_FOUND)
  set(CP2K_OPENBLAS64_ROOT ${CP2K_OPENBLAS_ROOT})
  set(CP2K_OPENBLA_THREADS_ROOT ${CP2K_OPENBLAS_ROOT})
  cp2k_find_libraries(OPENBLAS "openblas")
  cp2k_find_libraries(OPENBLAS64 "openblas64")
  cp2k_find_libraries(OPENBLAS_THREADS "openblas_threads;openblas_omp")
  cp2k_find_libraries(OPENBLAS_THREADS64 "openblas64_threads;openblas64_omp")
endif()

cp2k_include_dirs(OPENBLAS "cblas.h")

# check if found
if(CP2K_OPENBLAS_INCLUDE_DIRS)
  find_package_handle_standard_args(
    OpenBLAS REQUIRED_VARS CP2K_OPENBLAS_INCLUDE_DIRS
                           CP2K_OPENBLAS_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(OpenBLAS
                                    REQUIRED_VARS CP2K_OPENBLAS_LINK_LIBRARIES)
endif()

# add target to link against
if(CP2K_OPENBLAS_FOUND)
  if(NOT TARGET cp2k::BLAS::OpenBLAS::openblas)
    add_library(cp2k::BLAS::OpenBLAS::openblas INTERFACE IMPORTED)
  endif()
  set_property(
    TARGET cp2k::BLAS::OpenBLAS::openblas
    PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_OPENBLAS_LINK_LIBRARIES})
  if(CP2K_OPENBLAS_INCLUDE_DIRS)
    set_property(
      TARGET cp2k::BLAS::OpenBLAS::openblas
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CP2K_OPENBLAS_INCLUDE_DIRS})
  endif()
  if(NOT TARGET cp2k::BLAS::OpenBLAS::blas)
    add_library(cp2k::BLAS::OpenBLAS::blas ALIAS cp2k::BLAS::OpenBLAS::openblas)
  endif()
  set(CP2K_BLAS_VENDOR "OpenBLAS")
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_BLAS_VENDOR CP2K_OPENBLAS_FOUND
                 CP2K_OPENBLAS_LINK_LIBRARIES CP2K_OPENBLAS_INCLUDE_DIRS)
