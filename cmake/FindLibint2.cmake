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

cp2k_set_default_paths(LIBINT2 "Libint2")

find_package(PkgConfig REQUIRED)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBINT2 IMPORTED_TARGET GLOBAL libint2)
endif()

if(NOT CP2K_LIBINT2_FOUND)
  cp2k_find_libraries(LIBINT2 int2)
endif()

if(NOT CP2K_LIBINT2_INCLUDE_DIRS)
  cp2k_include_dirs(LIBINT2 "libint2.h;libint2/atom.h")
endif()

find_file(
  CP2K_LIBINT2_MOD_FILE
  NAMES "libint_f.mod"
  PATHS "${CP2K_LIBINT2_ROOT}/include" "${CP2K_LIBINT2_ROOT}/include/libint2")

if(NOT CP2K_LIBINT2_MOD_FILE)
  message(FATAL_ERROR "Libint2 : Fortran support is missing")
endif()

find_package_handle_standard_args(
  Libint2 CP2K_LIBINT2_FOUND CP2K_LIBINT2_INCLUDE_DIRS
  CP2K_LIBINT2_LINK_LIBRARIES)

if(CP2K_LIBINT2_FOUND)
  if(NOT TARGET cp2k::Libint2::int2)
    add_library(cp2k::Libint2::int2 INTERFACE IMPORTED)
  endif()

  if(CP2K_LIBINT2_INCLUDE_DIRS)
    set_target_properties(
      cp2k::Libint2::int2 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                     "${CP2K_LIBINT2_INCLUDE_DIRS}")
  endif()

  set_target_properties(
    cp2k::Libint2::int2 PROPERTIES INTERFACE_LINK_LIBRARIES
                                   ${CP2K_LIBINT2_LINK_LIBRARIES})
endif()

mark_as_advanced(CP2K_LIBINT2_FOUND CP2K_LIBINT2_LINK_LIBRARIES
                 CP2K_LIBINT2_INCLUDE_DIRS)
