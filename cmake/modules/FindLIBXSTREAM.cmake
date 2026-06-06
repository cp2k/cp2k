#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)
find_package(PkgConfig QUIET)

cp2k_set_default_paths(LIBXSTREAM "LIBXSTREAM")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBXSTREAM QUIET IMPORTED_TARGET GLOBAL libxstream)
  if(CP2K_LIBXSTREAM_PREFIX)
    set(CP2K_LIBXSTREAM_ROOT "${CP2K_LIBXSTREAM_PREFIX}")
  endif()
endif()

set(_cp2k_libxstream_hints
    "${CP2K_LIBXSTREAM_ROOT}" "${DBCSR_LIBXSTREAMROOT}" "$ENV{LIBXSTREAMROOT}"
    "${CMAKE_SOURCE_DIR}/../libxstream" "$ENV{HOME}/libxstream")

if(NOT CP2K_LIBXSTREAM_INCLUDE_DIRS)
  find_path(
    CP2K_LIBXSTREAM_INCLUDE_DIRS
    NAMES libxstream.h
    HINTS ${_cp2k_libxstream_hints}
    PATH_SUFFIXES include)
endif()

if(NOT CP2K_LIBXSTREAM_LINK_LIBRARIES)
  set(_libxstream_suffixes_save ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
  find_library(
    CP2K_LIBXSTREAM_LINK_LIBRARIES
    NAMES xstream
    HINTS ${_cp2k_libxstream_hints}
    PATH_SUFFIXES lib lib64)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_libxstream_suffixes_save})
endif()

find_file(
  CP2K_LIBXSTREAM_OPENCL_SCRIPT
  NAMES tool_opencl.sh
  HINTS ${_cp2k_libxstream_hints}
  PATH_SUFFIXES scripts)

find_file(
  CP2K_LIBXSTREAM_OPENCL_COMMON
  NAMES libxstream_atomics.h
  HINTS ${_cp2k_libxstream_hints} "${CP2K_LIBXSTREAM_INCLUDE_DIRS}"
  PATH_SUFFIXES include/opencl opencl)

foreach(_root IN LISTS _cp2k_libxstream_hints)
  if(NOT CP2K_LIBXSTREAM_SAMPLES_DIR AND EXISTS "${_root}/samples")
    set(CP2K_LIBXSTREAM_SAMPLES_DIR "${_root}/samples")
  endif()
endforeach()

find_package_handle_standard_args(
  LIBXSTREAM DEFAULT_MSG CP2K_LIBXSTREAM_INCLUDE_DIRS
  CP2K_LIBXSTREAM_LINK_LIBRARIES CP2K_LIBXSTREAM_OPENCL_SCRIPT
  CP2K_LIBXSTREAM_OPENCL_COMMON)

if(LIBXSTREAM_FOUND AND NOT TARGET cp2k::LIBXSTREAM)
  add_library(cp2k::LIBXSTREAM INTERFACE IMPORTED)
  if(TARGET PkgConfig::CP2K_LIBXSTREAM)
    target_link_libraries(cp2k::LIBXSTREAM INTERFACE PkgConfig::CP2K_LIBXSTREAM)
  else()
    set_target_properties(
      cp2k::LIBXSTREAM
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXSTREAM_INCLUDE_DIRS}"
                 INTERFACE_LINK_LIBRARIES "${CP2K_LIBXSTREAM_LINK_LIBRARIES}")
  endif()
endif()

mark_as_advanced(
  CP2K_LIBXSTREAM_INCLUDE_DIRS CP2K_LIBXSTREAM_LINK_LIBRARIES
  CP2K_LIBXSTREAM_OPENCL_SCRIPT CP2K_LIBXSTREAM_OPENCL_COMMON
  CP2K_LIBXSTREAM_SAMPLES_DIR)
