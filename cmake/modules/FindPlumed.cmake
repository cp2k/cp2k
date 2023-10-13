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
find_package(PkgConfig)

cp2k_set_default_paths(PLUMED "Plumed")

# First try with pkg
if(PKG_CONFIG_FOUND)
  # plumed has a pkg-config module
  pkg_search_module(CP2K_PLUMED IMPORTED_TARGET GLOBAL plumed plumedInternals)
endif()

if(NOT ${CP2K_PLUMED_FOUND})
  cp2k_find_libraries(PLUMED "plumed")
  cp2k_include_dirs(PLUMED "plumed.h plumed/plumed.h")
endif()

if(CP2K_PLUMED_INCLUDE_DIRS)
  find_package_handle_standard_args(Plumed DEFAULT_MSG CP2K_PLUMED_INCLUDE_DIRS
                                    CP2K_PLUMED_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Plumed DEFAULT_MSG
                                    CP2K_PLUMED_LINK_LIBRARIES)
endif()

if(CP2K_PLUMED_FOUND)
  if(NOT TARGET cp2k::plumed::plumed)
    add_library(cp2k::plumed::plumed INTERFACE IMPORTED)
  endif()
  set_target_properties(
    cp2k::plumed::plumed
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_PLUMED_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_PLUMED_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_PLUMED_LINK_LIBRARIES CP2K_PLUMED_INCLUDE_DIRS
                 CP2K_PLUMED_FOUND)
