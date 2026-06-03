#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)

# Probe: user override > DBCSR hint > sibling of LIBXS > environment > common
# paths
if(NOT LIBXSMMROOT)
  foreach(_dir
          "${DBCSR_LIBXSMMROOT}" "$ENV{LIBXSMMROOT}" "${LIBXSROOT}/../libxsmm"
          "${CMAKE_SOURCE_DIR}/../libxsmm" "$ENV{HOME}/libxsmm")
    if(EXISTS "${_dir}/include/libxsmm.h")
      set(LIBXSMMROOT "${_dir}")
      break()
    endif()
  endforeach()
endif()

if(LIBXSMMROOT)
  find_path(
    CP2K_LIBXSMM_INCLUDE_DIR libxsmm.h
    PATHS "${LIBXSMMROOT}/include"
    NO_DEFAULT_PATH)
  set(_libxsmm_suffixes_save ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
  find_library(
    CP2K_LIBXSMM_LIBRARY xsmm
    PATHS "${LIBXSMMROOT}/lib"
    NO_DEFAULT_PATH)
  if(NOT CP2K_LIBXSMM_LIBRARY)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_libxsmm_suffixes_save})
    find_library(
      CP2K_LIBXSMM_LIBRARY xsmm
      PATHS "${LIBXSMMROOT}/lib"
      NO_DEFAULT_PATH)
  endif()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_libxsmm_suffixes_save})
endif()

find_package_handle_standard_args(LIBXSMM DEFAULT_MSG CP2K_LIBXSMM_INCLUDE_DIR
                                  CP2K_LIBXSMM_LIBRARY)

if(LIBXSMM_FOUND AND NOT TARGET cp2k::LIBXSMM)
  add_library(cp2k::LIBXSMM UNKNOWN IMPORTED)
  set_target_properties(
    cp2k::LIBXSMM
    PROPERTIES IMPORTED_LOCATION "${CP2K_LIBXSMM_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXSMM_INCLUDE_DIR}")
endif()

mark_as_advanced(CP2K_LIBXSMM_INCLUDE_DIR CP2K_LIBXSMM_LIBRARY)
