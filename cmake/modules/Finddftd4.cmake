#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(DFTD4 "dftd4")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(DFTD4 QUIET IMPORTED_TARGET GLOBAL dftd4)
endif()

if(NOT DFTD4_FOUND)
  cp2k_find_libraries(DFTD4 dftd4)
endif()

if(NOT DFTD4_INCLUDE_DIRS)
  cp2k_include_dirs(DFTD4 "dftd4.h")
endif()

if(DFTD4_INCLUDE_DIRS)
  find_package_handle_standard_args(dftd4 DEFAULT_MSG DFTD4_INCLUDE_DIRS
                                    DFTD4_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(dftd4 DEFAULT_MSG DFTD4_LINK_LIBRARIES)
endif()

if(DFTD4_FOUND)
  if(NOT TARGET dftd4::dftd4)
    add_library(dftd4::dftd4 INTERFACE IMPORTED)
  endif()
  set_target_properties(
    dftd4::dftd4
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${DFTD4_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${DFTD4_LINK_LIBRARIES}")
endif()

mark_as_advanced(DFTD4_LINK_LIBRARIES DFTD4_INCLUDE_DIRS DFTD4_FOUND)
