#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(TBLITE "tblite")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(TBLITE QUIET IMPORTED_TARGET GLOBAL tblite)
endif()

if(NOT TBLITE_FOUND)
  cp2k_find_libraries(TBLITE tblite)
endif()

if(NOT TBLITE_INCLUDE_DIRS)
  cp2k_include_dirs(TBLITE "tblite.h")
endif()

if(TBLITE_INCLUDE_DIRS)
  find_package_handle_standard_args(tblite DEFAULT_MSG TBLITE_INCLUDE_DIRS
                                    TBLITE_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(tblite DEFAULT_MSG TBLITE_LINK_LIBRARIES)
endif()

if(TBLITE_FOUND)
  if(NOT TARGET tblite::tblite)
    add_library(tblite::tblite INTERFACE IMPORTED)
  endif()
  set_target_properties(
    tblite::tblite
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TBLITE_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${TBLITE_LINK_LIBRARIES}")
endif()

mark_as_advanced(TBLITE_LINK_LIBRARIES TBLITE_INCLUDE_DIRS TBLITE_FOUND)
