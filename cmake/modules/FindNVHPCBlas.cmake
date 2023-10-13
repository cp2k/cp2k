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

cp2k_set_default_paths(CP2K_NVHPC "NVHPC")

find_library(
  CP2K_NVHPC_BLAS_LP64
  NAMES blas_lp64
  PATHS "${CP2K_NVHPC_ROOT}"
  PATH_SUFFIXES "lib" "lib64")

find_library(
  CP2K_NVHPC_BLAS_ILP64
  NAMES blas_ilp64
  PATHS "${CP2K_NVHPC_ROOT}"
  PATH_SUFFIXES "lib" "lib64")

find_library(
  CP2K_NVHPC_LAPACK_LP64
  NAMES lapack_lp64
  PATHS "${CP2K_NVHPC_ROOT}"
  PATH_SUFFIXES "lib" "lib64")

find_library(
  CP2K_NVHPC_LAPACK_ILP64
  NAMES lapack_ilp64
  PATHS "${CP2K_NVHPC_ROOT}"
  PATH_SUFFIXES "lib" "lib64")

find_path(
  CP2K_NVHPC_BLAS_INCLUDE_DIRS_lp64
  NAMES cblas.h
  PATHS "${CP2K_NVHPC_ROOT}"
  HINTS "${CP2K_NVHPC_ROOT}"
  PATH_SUFFIXES "include" "include/lp64" "lp64")

find_path(
  CP2K_NVHPC_BLAS_INCLUDE_DIRS_ilp64
  NAMES cblas.h
  PATHS "${CP2K_NVHPC_ROOT}"
  HINTS "${CP2K_NVHPC_ROOT}"
  PATH_SUFFIXES "include" "include/ilp64" "ilp64")

find_package_handle_standard_args(
  NVHPCBlas
  DEFAULT_MSG
  CP2K_NVHPC_INCLUDE_DIRS_ipl64
  CP2K_NVHPC_BLAS_INCLUDE_DIRS_lp64
  CP2K_NVHPC_BLAS_ILP64
  CP2K_NVHPC_BLAS_LP64
  CP2K_NVHPC_LAPACK_ILP64
  CP2K_NVHPC_LAPACK_LP64)

set(CP2K_BLAS_VENDOR "NVHPCBlas")
set(CP2K_NVHPCBLAS_FOUND "ON")

if(NOT TARGET cp2k::BLAS::NVHPCBlas::nvhpcblas)
  add_library(cp2k::BLAS::NVHPCBlas::nvhpcblas INTERFACE IMPORTED)
  add_library(cp2k::BLAS::NVHPCBlas::blas ALIAS
              cp2k::BLAS::NVHPCBlas::nvhpcblas)
endif()

if(CP2K_BLAS_INTERFACE MATCHES "64bits")
  set(CP2K_NVHPC_BLAS_LINK_LIBRARIES
      "${CP2K_NVHPC_LAPACK_ILP64} ${CP2K_NVHPC_BLAS_ILP64}")
  set(CP2K_NVHPC_BLAS_INCLUDE_DIRS "${CP2K_NVHPC_INCLUDE_DIRS_ipl64}")
else()
  set(CP2K_NVHPC_BLAS_LINK_LIBRARIES "${CP2K_NVHPC_LAPACK_LP64}
        ${CP2K_NVHPC_BLAS_LP64}")
  set(CP2K_NVHPC_BLAS_INCLUDE_DIRS "${CP2K_NVHPC_INCLUDE_DIRS_pl64}")
endif()

set_target_properties(
  cp2k::BLAS::NVHPCBlas::nvhpcblas
  PROPERTIES INTERFACE_LINK_LIBRARIES "${CP2K_NVHPC_BLAS_LINK_LIBRARIES}"
             INTERFACE_INCLUDE_DIRECTORIES "${CP2K_NVHPC_BLAS_INCLUDE_DIRS}")

mark_as_advanced(CP2K_NVHPCBLAS_FOUND CP2K_NVHPC_BLAS_INCLUDE_DIRS
                 CP2K_NVHPC_BLAS_LINK_LIBRARIES CP2K_BLAS_VENDOR)
