#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

if(NOT POLICY CMP0074)
  set(_GenericBLAS_PATHS ${GenericBLAS_ROOT} $ENV{GenericBLAS_ROOT})
endif()

find_library(
  GenericBLAS_LIBRARIES
  NAMES "blas"
  HINTS ${_GenericBLAS_PATHS})
find_library(
  # optinally look for cblas library - not required
  GenericBLAS_CBLAS_LIBRARIES
  NAMES "cblas"
  HINTS ${_GenericBLAS_PATHS})
find_path(
  GenericBLAS_INCLUDE_DIRS
  NAMES "cblas.h"
  HINTS ${_GenericBLAS_PATHS})

# check if found
include(FindPackageHandleStandardArgs)
if(GenericBLAS_INCLUDE_DIRS)
  find_package_handle_standard_args(
    GenericBLAS REQUIRED_VARS GenericBLAS_INCLUDE_DIRS GenericBLAS_LIBRARIES)
else()
  find_package_handle_standard_args(GenericBLAS
                                    REQUIRED_VARS GenericBLAS_LIBRARIES)
endif()

if(GenericBLAS_CBLAS_LIBRARIES)
  list(APPEND GenericBLAS_LIBRARIES ${GenericBLAS_CBLAS_LIBRARIES})
endif()

# add target to link against
if(GenericBLAS_FOUND)

  if(NOT TARGET cp2k::BLAS::GenericBLAS::blas)
    add_library(cp2k::BLAS::GenericBLAS::blas INTERFACE IMPORTED)
  endif()
  set_property(TARGET cp2k::BLAS::GenericBLAS::blas
               PROPERTY INTERFACE_LINK_LIBRARIES ${GenericBLAS_LIBRARIES})
  set_property(
    TARGET cp2k::BLAS::GenericBLAS::blas PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                  ${GenericBLAS_INCLUDE_DIRS})
endif()

# prevent clutter in cache
mark_as_advanced(GenericBLAS_FOUND GenericBLAS_LIBRARIES
                 GenericBLAS_INCLUDE_DIRS GenericBLAS_CBLAS_LIBRARIES)
