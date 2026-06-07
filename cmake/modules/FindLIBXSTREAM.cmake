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

cp2k_set_default_paths(LIBXSTREAM "LIBXSTREAM")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBXSTREAM QUIET IMPORTED_TARGET GLOBAL libxstream)
endif()

set(_cp2k_libxstream_hints)
if(CP2K_LIBXSTREAM_ROOT)
  list(APPEND _cp2k_libxstream_hints "${CP2K_LIBXSTREAM_ROOT}")
endif()
if(DEFINED ENV{LIBXSTREAMROOT})
  list(APPEND _cp2k_libxstream_hints "$ENV{LIBXSTREAMROOT}")
endif()

if(NOT CP2K_LIBXSTREAM_INCLUDE_DIR)
  find_path(
    CP2K_LIBXSTREAM_INCLUDE_DIR
    NAMES libxstream.h
    HINTS ${_cp2k_libxstream_hints}
    PATH_SUFFIXES include)
endif()

if(NOT CP2K_LIBXSTREAM_LIBRARY)
  find_library(
    CP2K_LIBXSTREAM_LIBRARY
    NAMES xstream
    HINTS ${_cp2k_libxstream_hints}
    PATH_SUFFIXES lib lib64)
endif()

if(NOT CP2K_LIBXSTREAM_OPENCL_SCRIPT)
  find_file(
    CP2K_LIBXSTREAM_OPENCL_SCRIPT
    NAMES tool_opencl.sh
    HINTS ${_cp2k_libxstream_hints}
    PATH_SUFFIXES scripts)
endif()

if(NOT CP2K_LIBXSTREAM_OPENCL_COMMON)
  find_file(
    CP2K_LIBXSTREAM_OPENCL_COMMON
    NAMES libxstream_atomics.h
    HINTS ${_cp2k_libxstream_hints}
    PATH_SUFFIXES include/opencl opencl)
endif()

if(NOT CP2K_LIBXSTREAM_SAMPLES_DIR)
  foreach(_root IN LISTS _cp2k_libxstream_hints)
    if(EXISTS "${_root}/samples")
      set(CP2K_LIBXSTREAM_SAMPLES_DIR "${_root}/samples")
      break()
    endif()
  endforeach()
endif()

# Reconcile pkg-config variables with local naming.
if(NOT CP2K_LIBXSTREAM_INCLUDE_DIR AND CP2K_LIBXSTREAM_INCLUDE_DIRS)
  list(GET CP2K_LIBXSTREAM_INCLUDE_DIRS 0 CP2K_LIBXSTREAM_INCLUDE_DIR)
endif()

if(NOT CP2K_LIBXSTREAM_LIBRARY AND CP2K_LIBXSTREAM_LINK_LIBRARIES)
  list(GET CP2K_LIBXSTREAM_LINK_LIBRARIES 0 CP2K_LIBXSTREAM_LIBRARY)
endif()

if(NOT CP2K_LIBXSTREAM_ROOT AND CP2K_LIBXSTREAM_PREFIX)
  set(CP2K_LIBXSTREAM_ROOT "${CP2K_LIBXSTREAM_PREFIX}")
endif()

find_package_handle_standard_args(
  LIBXSTREAM
  REQUIRED_VARS CP2K_LIBXSTREAM_INCLUDE_DIR CP2K_LIBXSTREAM_LIBRARY
                CP2K_LIBXSTREAM_OPENCL_SCRIPT CP2K_LIBXSTREAM_OPENCL_COMMON)

if(LIBXSTREAM_FOUND AND NOT TARGET cp2k::libxstream)
  add_library(cp2k::libxstream INTERFACE IMPORTED)

  if(TARGET PkgConfig::CP2K_LIBXSTREAM)
    target_link_libraries(cp2k::libxstream INTERFACE PkgConfig::CP2K_LIBXSTREAM)
  else()
    set_target_properties(
      cp2k::libxstream
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXSTREAM_INCLUDE_DIR}"
                 INTERFACE_LINK_LIBRARIES "${CP2K_LIBXSTREAM_LIBRARY}")
  endif()
endif()

mark_as_advanced(
  CP2K_LIBXSTREAM_INCLUDE_DIR CP2K_LIBXSTREAM_LIBRARY
  CP2K_LIBXSTREAM_OPENCL_SCRIPT CP2K_LIBXSTREAM_OPENCL_COMMON
  CP2K_LIBXSTREAM_SAMPLES_DIR)
