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

cp2k_set_default_paths(ATLAS "Atlas")

cp2k_find_libraries(ATLAS "atlas")
cp2k_include_dirs(ATLAS "cblas.h atlas/cblas.h")
# check if found
find_package_handle_standard_args(Atlas REQUIRED_VARS CP2K_ATLAS_INCLUDE_DIRS
                                                      CP2K_ATLAS_LINK_LIBRARIES)

# add target to link against
if(CP2K_ATLAS_FOUND)
  if(NOT TARGET cp2k::BLAS::ATLAS::atlas)
    add_library(cp2k::BLAS::ATLAS::atlas INTERFACE IMPORTED)
  endif()
  set_property(TARGET cp2k::BLAS::ATLAS::atlas
               PROPERTY INTERFACE_LINK_LIBRARIES ${CP2K_ATLAS_LINK_LIBRARIES})
  set_property(
    TARGET cp2k::BLAS::ATLAS::atlas PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                             ${CP2K_ATLAS_INCLUDE_DIRS})
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_ATLAS_FOUND CP2K_ATLAS_LINK_LIBRARIES
                 CP2K_ATLAS_INCLUDE_DIRS)
