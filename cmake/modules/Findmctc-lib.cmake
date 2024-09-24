#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(MCTC-LIB "mctc-lib")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(MCTC-LIB QUIET IMPORTED_TARGET GLOBAL mctc-lib)
endif()

if(NOT MCTC-LIB_FOUND)
  cp2k_find_libraries(MCTC-LIB mctc-lib)
endif()

if(MCTC-LIB_INCLUDE_DIRS)
  find_package_handle_standard_args(mctc-lib DEFAULT_MSG MCTC-LIB_INCLUDE_DIRS
                                    MCTC-LIB_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(mctc-lib DEFAULT_MSG
                                    MCTC-LIB_LINK_LIBRARIES)
endif()

if(MCTC-LIB_FOUND)
  if(NOT TARGET mctc-lib::mctc-lib)
    add_library(mctc-lib::mctc-lib INTERFACE IMPORTED)
  endif()
  set_target_properties(
    mctc-lib::mctc-lib
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MCTC-LIB_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${MCTC-LIB_LINK_LIBRARIES}")
endif()

mark_as_advanced(MCTC-LIB_LINK_LIBRARIES MCTC-LIB_INCLUDE_DIRS MCTC-LIB_FOUND)
