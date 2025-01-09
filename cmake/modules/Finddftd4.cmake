#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)
find_package(PkgConfig)

cp2k_set_default_paths(DFTD4 "dftd4")

if(PKG_CONFIG_FOUND)
  pkg_search_module(CP2K_DFTD4 IMPORTED_TARGET GLOBAL dftd4)
endif()

if(NOT ${CP2K_DFTD4_FOUND})
  cp2k_find_libraries(DFTD4 "dftd4")
  cp2k_include_dirs(DFTD4 "dftd4.h")
endif()

if(CP2K_DFTD4_INCLUDE_DIRS)
  find_package_handle_standard_args(dftd4 DEFAULT_MSG CP2K_DFTD4_INCLUDE_DIRS
                                    CP2K_DFTD4_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(dftd4 DEFAULT_MSG CP2K_DFTD4_LINK_LIBRARIES)
endif()

if(CP2K_DFTD4_FOUND)
  if(NOT TARGET cp2k::dftd4)
    add_library(cp2k::dftd4 INTERFACE IMPORTED)
  endif()
  set_target_properties(
    cp2k::dftd4
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_DFTD4_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_DFTD4_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_DFTD4_LINK_LIBRARIES CP2K_DFTD4_INCLUDE_DIRS
                 CP2K_DFTD4_FOUND)
