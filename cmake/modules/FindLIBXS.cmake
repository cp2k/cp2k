#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
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
  if(CP2K_LIBXS_PREFIX)
    set(CP2K_LIBXS_ROOT "${CP2K_LIBXS_PREFIX}")
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

set(_cp2k_libxs_required_headers libxs.h libxs_gemm.h libxs_malloc.h
                                 libxs_timer.h)
set(_cp2k_libxs_headers_found TRUE)
foreach(_header IN LISTS _cp2k_libxs_required_headers)
  if(NOT CP2K_LIBXS_INCLUDE_DIRS OR NOT EXISTS
                                    "${CP2K_LIBXS_INCLUDE_DIRS}/${_header}")
    set(_cp2k_libxs_headers_found FALSE)
  endif()
endforeach()

if(NOT _cp2k_libxs_headers_found)
  set(CP2K_LIBXS_INCLUDE_DIRS "CP2K_LIBXS_INCLUDE_DIRS-NOTFOUND")
endif()

if(NOT CP2K_LIBXS_LINK_LIBRARIES)
  set(_libxs_suffixes_save ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
  find_library(
    CP2K_LIBXS_LINK_LIBRARIES
    NAMES xs
    HINTS "${CP2K_LIBXS_ROOT}" "${DBCSR_LIBXSROOT}"
          "${CMAKE_SOURCE_DIR}/../libxs" "$ENV{HOME}/libxs" "/opt/libxs"
    PATH_SUFFIXES lib lib64)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_libxs_suffixes_save})
endif()

find_package_handle_standard_args(LIBXS DEFAULT_MSG CP2K_LIBXS_INCLUDE_DIRS
                                  CP2K_LIBXS_LINK_LIBRARIES)

if(LIBXS_FOUND AND NOT TARGET cp2k::LIBXS)
  add_library(cp2k::LIBXS INTERFACE IMPORTED)
  if(TARGET PkgConfig::CP2K_LIBXS)
    target_link_libraries(cp2k::LIBXS INTERFACE PkgConfig::CP2K_LIBXS)
  else()
    set_target_properties(
      cp2k::LIBXS
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXS_INCLUDE_DIRS}"
                 INTERFACE_LINK_LIBRARIES "${CP2K_LIBXS_LINK_LIBRARIES}")
  endif()
endif()

mark_as_advanced(CP2K_LIBXS_INCLUDE_DIRS CP2K_LIBXS_LINK_LIBRARIES)
