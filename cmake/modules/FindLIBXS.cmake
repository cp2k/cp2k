#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(CheckIncludeFiles)
include(cp2k_utils)
find_package(PkgConfig QUIET)

cp2k_set_default_paths(LIBXS "LIBXS")

# Prefer pkg-config when available.  This covers package-manager installs and
# keeps transitive flags in one place.
if(PKG_CONFIG_FOUND)
  if(BUILD_SHARED_LIBS)
    pkg_check_modules(CP2K_LIBXS QUIET IMPORTED_TARGET GLOBAL libxs-shared)
  else()
    pkg_check_modules(CP2K_LIBXS QUIET IMPORTED_TARGET GLOBAL libxs-static)
  endif()
  if(NOT CP2K_LIBXS_FOUND)
    pkg_check_modules(CP2K_LIBXS QUIET IMPORTED_TARGET GLOBAL libxs)
  endif()
endif()

# Fallback to regular CMake search paths plus explicit hints.  Do not restrict
# the search with NO_DEFAULT_PATH, otherwise /usr, /usr/local, Spack prefixes,
# and CMAKE_PREFIX_PATH are ignored.
if(NOT CP2K_LIBXS_INCLUDE_DIRS)
  find_path(
    CP2K_LIBXS_INCLUDE_DIRS
    NAMES libxs.h
    HINTS "${CP2K_LIBXS_ROOT}" "${DBCSR_LIBXSROOT}"
          "${CMAKE_SOURCE_DIR}/../libxs" "$ENV{HOME}/libxs" "/opt/libxs"
    PATH_SUFFIXES include)
endif()

# Multiple headers in one shot. Sanity check of their existence
check_include_files(
  "libxs/libxs_gemm.h;libxs/libxs_malloc.h;libxs/libxs_timer.h"
  HAVE_LIBXS_HEADERS)

if(NOT CP2K_LIBXS_LINK_LIBRARIES)
  # Assume CP2K_LIBXS_ROOT is an optional explicit prefix; otherwise rely on
  # CMAKE_PREFIX_PATH.
  set(_libxs_hints)
  if(CP2K_LIBXS_ROOT)
    list(APPEND _libxs_hints "${CP2K_LIBXS_ROOT}")
  endif()

  find_library(
    CP2K_LIBXS_LIBRARY
    NAMES xs
    HINTS ${_libxs_hints}
    PATH_SUFFIXES lib lib64)
endif()

# check that all requirements are met
find_package_handle_standard_args(LIBXS DEFAULT_MSG CP2K_LIBXS_INCLUDE_DIRS
                                  CP2K_LIBXS_LINK_LIBRARIES HAVE_LIBXS_HEADERS)

if(LIBXS_FOUND AND NOT TARGET cp2k::libxs)
  if(TARGET PkgConfig::CP2K_LIBXS)
    add_library(cp2k::libxs ALIAS PkgConfig::CP2K_LIBXS)
  else()
    add_library(cp2k::libxs INTERFACE IMPORTED)
    set_target_properties(
      cp2k::libxs
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXS_INCLUDE_DIRS}"
                 INTERFACE_LINK_LIBRARIES "${CP2K_LIBXS_LINK_LIBRARIES}")
  endif()
endif()

mark_as_advanced(CP2K_LIBXS_INCLUDE_DIRS CP2K_LIBXS_LINK_LIBRARIES)
