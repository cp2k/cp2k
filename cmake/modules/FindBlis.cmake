#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

find_package(PkgConfig)
include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(BLIS "BLIS")

if(DEFINED AOCL_ROOT)
  list(CP2K_BLIS_ROOT "${AOCL_ROOT}" "$ENV{AOCL_ROOT}")
endif()

# one day blis will have a pkg-config file
if(PKG_CONFIG_FOUND)
  pkg_check_modules(BLIS IMPORTED_TARGET GLOBAL blis)
endif()

if(NOT CP2K_BLIS_FOUND)
  cp2k_find_libraries(BLIS "blis")
endif()

if(NOT CP2K_BLIS_INCLUDE_DIRS)
  cp2k_include_dirs(BLIS "blis.h")
endif()

# check if found
if(CP2K_BLIS_INCLUDE_DIRS)
  find_package_handle_standard_args(
    BLIS REQUIRED_VARS CP2K_BLIS_FOUND CP2K_BLIS_INCLUDE_DIRS
                       CP2K_BLIS_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Blis REQUIRED_VARS CP2K_BLIS_FOUND
                                                       CP2K_BLIS_LINK_LIBRARIES)
endif()

# add target to link against
if(CP2K_BLIS_FOUND)
  if(NOT TARGET cp2k::BLAS::Blis::blis)
    add_library(cp2k::BLAS::Blis::blis INTERFACE IMPORTED)
    add_library(cp2k::BLAS::Blis::blas alias cp2k::BLAS::Blis::blis)
  endif()

  set_property(TARGET cp2k::BLAS::Blis::blis
               PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_BLIS_LINK_LIBRARIES})

  if(BLIS_INCLUDE_DIRS)
    set_property(
      TARGET cp2k::BLAS::Blis::blis PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                             ${CP2K_BLIS_INCLUDE_DIRS})
  endif()
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_BLIS_FOUND CP2K_BLIS_LINK_LIBRARIES
                 CP2K_BLIS_INCLUDE_DIRS)
