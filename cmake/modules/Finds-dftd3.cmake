#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(S-DFTD3 "s-dftd3")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(S-DFTD3 QUIET IMPORTED_TARGET GLOBAL s-dftd3)
endif()

if(NOT S-DFTD3_FOUND)
  cp2k_find_libraries(S-DFTD3 s-dftd3)
endif()

if(NOT S-DFTD3_INCLUDE_DIRS)
  cp2k_include_dirs(S-DFTD3 "s-dftd3.h")
endif()

if(S-DFTD3_INCLUDE_DIRS)
  find_package_handle_standard_args(s-dftd3 DEFAULT_MSG S-DFTD3_INCLUDE_DIRS
                                    S-DFTD3_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(s-dftd3 DEFAULT_MSG S-DFTD3_LINK_LIBRARIES)
endif()

if(S-DFTD3_FOUND)
  if(NOT TARGET s-dftd3::s-dftd3)
    add_library(s-dftd3::s-dftd3 INTERFACE IMPORTED GLOBAL)
  endif()
  get_filename_component(S-DFTD3_LINK_LIB "${S-DFTD3_LINK_LIBRARIES}" PATH)
  set_target_properties(
    s-dftd3::s-dftd3
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${S-DFTD3_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${S-DFTD3_LINK_LIBRARIES}")
endif()

mark_as_advanced(S-DFTD3_LINK_LIB S-DFTD3_LINK_LIBRARIES S-DFTD3_INCLUDE_DIRS
                 S-DFTD3_FOUND)
