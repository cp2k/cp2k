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
  pkg_search_module(PKGC_DFTD4 IMPORTED_TARGET GLOBAL dftd4)
  pkg_search_module(PKGC_MCTC mctc-lib)
  pkg_search_module(PKGC_MULTICHARGE multicharge)
endif()

if(dftd4_FIND_REQUIRED AND NOT PKGC_DFTD4_FOUND)
  message(FATAL_ERROR "Unable to find DFTD4 and/or its dependencies.")
endif()

set(CP2K_DFTD4_FOUND TRUE)

find_path(
  CP2K_DFTD4_INCLUDE_DIR
  NAMES dftd4.h
  PATHS ${PKGC_DFTD4_INCLUDE_DIRS})

find_path(
  CP2K_DFTD4_MOD_DIR
  NAMES dftd4.mod
  PATHS ${PKGC_DFTD4_INCLUDE_DIRS})
message(STATUS "DFTD4 mod dir: ${CP2K_DFTD4_MOD_DIR}")

find_library(
  CP2K_DFTD4_LINK_LIBRARIES
  NAMES ${PKGC_DFTD4_LIBRARIES}
  PATHS ${PKGC_DFTD4_LIBRARY_DIRS}
  DOC "dftd4 libraries")
find_library(
  CP2K_MCTC_LINK_LIBRARIES
  NAMES ${PKGC_MCTC_LIBRARIES}
  PATHS ${PKGC_MCTC_LIBRARY_DIRS}
  DOC "mctc libraries")
find_library(
  CP2K_MULTICHARGE_LINK_LIBRARIES
  NAMES ${PKGC_MULTICHARGE_LIBRARIES}
  PATHS ${PKGC_MULTICHARGE_LIBRARY_DIRS}
  DOC "multicharge libraries")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(dftd4 DEFAULT_MSG CP2K_DFTD4_INCLUDE_DIR
                                  CP2K_DFTD4_LINK_LIBRARIES)

if(CP2K_DFTD4_FOUND)
  add_library(cp2k::dftd4 INTERFACE IMPORTED)
  set_target_properties(
    cp2k::dftd4
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES
      "${CP2K_DFTD4_INCLUDE_DIR};${CP2K_DFTD4_MOD_DIR}"
      INTERFACE_LINK_LIBRARIES
      "${CP2K_DFTD4_LINK_LIBRARIES};${CP2K_MCTC_LINK_LIBRARIES};${CP2K_MULTICHARGE_LINK_LIBRARIES}"
  )
endif()
