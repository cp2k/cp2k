#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)
find_package(PkgConfig QUIET)

cp2k_set_default_paths(LIBXSMM "LIBXSMM")

# Prefer pkg-config when available.
if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBXSMM QUIET IMPORTED_TARGET GLOBAL libxsmm)
endif()

# Optional explicit prefix hint.
set(_cp2k_libxsmm_hints)
if(CP2K_LIBXSMM_ROOT)
  list(APPEND _cp2k_libxsmm_hints "${CP2K_LIBXSMM_ROOT}")
endif()

# Fallback discovery when pkg-config did not populate usable values.
if(NOT CP2K_LIBXSMM_INCLUDE_DIR)
  find_path(
    CP2K_LIBXSMM_INCLUDE_DIR
    NAMES libxsmm.h
    HINTS ${_cp2k_libxsmm_hints}
    PATH_SUFFIXES include)
endif()

if(NOT CP2K_LIBXSMM_LIBRARY)
  find_library(
    CP2K_LIBXSMM_LIBRARY
    NAMES xsmm
    HINTS ${_cp2k_libxsmm_hints}
    PATH_SUFFIXES lib lib64)
endif()

# Reconcile pkg-config variables with the fallback naming used below.
if(NOT CP2K_LIBXSMM_INCLUDE_DIR AND CP2K_LIBXSMM_INCLUDE_DIRS)
  list(GET CP2K_LIBXSMM_INCLUDE_DIRS 0 CP2K_LIBXSMM_INCLUDE_DIR)
endif()

if(NOT CP2K_LIBXSMM_LIBRARY AND CP2K_LIBXSMM_LINK_LIBRARIES)
  list(GET CP2K_LIBXSMM_LINK_LIBRARIES 0 CP2K_LIBXSMM_LIBRARY)
endif()

find_package_handle_standard_args(LIBXSMM REQUIRED_VARS CP2K_LIBXSMM_INCLUDE_DIR
                                                        CP2K_LIBXSMM_LIBRARY)

if(LIBXSMM_FOUND AND NOT TARGET cp2k::libxsmm)
  add_library(cp2k::libxsmm INTERFACE IMPORTED)

  if(TARGET PkgConfig::CP2K_LIBXSMM)
    target_link_libraries(cp2k::libxsmm INTERFACE PkgConfig::CP2K_LIBXSMM)
  else()
    set_target_properties(
      cp2k::libxsmm
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXSMM_INCLUDE_DIR}"
                 INTERFACE_LINK_LIBRARIES "${CP2K_LIBXSMM_LIBRARY}")
  endif()
endif()

mark_as_advanced(CP2K_LIBXSMM_INCLUDE_DIR CP2K_LIBXSMM_LIBRARY)
