#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(TOML-F "toml-f")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(TOML-F QUIET IMPORTED_TARGET GLOBAL toml-f)
endif()

if(NOT TOML-F_FOUND)
  cp2k_find_libraries(TOML-F toml-f)
endif()

if(TOML-F_INCLUDE_DIRS)
  find_package_handle_standard_args(toml-f DEFAULT_MSG TOML-F_INCLUDE_DIRS
                                    TOML-F_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(toml-f DEFAULT_MSG TOML-F_LINK_LIBRARIES)
endif()

if(TOML-F_FOUND)
  if(NOT TARGET toml-f::toml-f)
    add_library(toml-f::toml-f INTERFACE IMPORTED GLOBAL)
  endif()
  get_filename_component(TOML-F_LINK_LIB "${TOML-F_LINK_LIBRARIES}" PATH)
  set_target_properties(
    toml-f::toml-f
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TOML-F_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${TOML-F_LINK_LIBRARIES}")
endif()

mark_as_advanced(TOML-F_LINK_LIB TOML-F_LINK_LIBRARIES TOML-F_INCLUDE_DIRS
                 TOML-F_FOUND)
