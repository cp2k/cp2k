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

cp2k_set_default_paths(LIBSPG "LibSPG")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBSPG IMPORTED_TARGET spglib)
endif()

if(NOT CP2K_LIBSPG_FOUND)
  cp2k_find_libraries(LIBSPG "symspg")
endif()

if(NOT DEFINED CP2K_LIBSPG_INCLUDE_DIRS)
  cp2k_include_dirs(LIBSPG "spglib.h")
endif()

if(CP2K_LIBSPG_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibSPG DEFAULT_MSG CP2K_LIBSPG_LINK_LIBRARIES CP2K_LIBSPG_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibSPG DEFAULT_MSG
                                    CP2K_LIBSPG_LINK_LIBRARIES)
endif()

if(CP2K_LIBSPG_FOUND AND NOT TARGET cp2k::LIBSPG::libspg)
  add_library(cp2k::LIBSPG::libspg INTERFACE IMPORTED)
  set_target_properties(
    cp2k::LIBSPG::libspg PROPERTIES INTERFACE_LINK_LIBRARIES
                                    "${CP2K_LIBSPG_LINK_LIBRARIES}")
  if(CP2K_LIBSPG_INCLUDE_DIRS)
    set_target_properties(
      cp2k::LIBSPG::libspg PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                      "${CP2K_LIBSPG_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBSPG_LINK_LIBRARIES)
mark_as_advanced(CP2K_LIBSPG_INCLUDE_DIRS)
mark_as_advanced(CP2K_LIBSPG_FOUND)
