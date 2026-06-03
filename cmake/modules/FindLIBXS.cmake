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

# Path probing: DBCSR hint > sibling of source > $HOME > /opt
if(CP2K_LIBXS_ROOT STREQUAL "/usr")
  if(DBCSR_LIBXSROOT AND EXISTS "${DBCSR_LIBXSROOT}/include/libxs.h")
    set(CP2K_LIBXS_ROOT "${DBCSR_LIBXSROOT}")
  else()
    foreach(_dir "${CMAKE_SOURCE_DIR}/../libxs" "$ENV{HOME}/libxs" "/opt/libxs")
      if(EXISTS "${_dir}/include/libxs.h")
        set(CP2K_LIBXS_ROOT "${_dir}")
        break()
      endif()
    endforeach()
  endif()
endif()

if(PKG_CONFIG_FOUND)
  if(BUILD_SHARED_LIBS)
    pkg_check_modules(CP2K_LIBXS IMPORTED_TARGET GLOBAL libxs-shared)
  else()
    pkg_check_modules(CP2K_LIBXS IMPORTED_TARGET GLOBAL libxs-static)
  endif()
  if(NOT CP2K_LIBXS_FOUND)
    pkg_check_modules(CP2K_LIBXS IMPORTED_TARGET GLOBAL libxs)
  endif()
endif()

if(NOT CP2K_LIBXS_FOUND)
  set(_libxs_suffixes_save ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib)
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
  cp2k_find_libraries(LIBXS xs)
  if(NOT CP2K_LIBXS_FOUND)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_libxs_suffixes_save})
    cp2k_find_libraries(LIBXS xs)
  endif()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_libxs_suffixes_save})
endif()

if(NOT CP2K_LIBXS_INCLUDE_DIRS)
  cp2k_include_dirs(LIBXS "libxs.h;include/libxs.h")
endif()

if(CP2K_LIBXS_INCLUDE_DIRS)
  find_package_handle_standard_args(LIBXS DEFAULT_MSG CP2K_LIBXS_INCLUDE_DIRS
                                    CP2K_LIBXS_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(LIBXS DEFAULT_MSG CP2K_LIBXS_LINK_LIBRARIES)
endif()

if(NOT TARGET cp2k::LIBXS::libxs)
  add_library(cp2k::LIBXS::libxs INTERFACE IMPORTED)
  if(CP2K_LIBXS_FOUND)
    if(CP2K_LIBXS_LIBRARY_DIRS)
      target_link_directories(cp2k::LIBXS::libxs INTERFACE
                              ${CP2K_LIBXS_LIBRARY_DIRS})
    endif()
    set_target_properties(
      cp2k::LIBXS::libxs PROPERTIES INTERFACE_LINK_LIBRARIES
                                    "${CP2K_LIBXS_LINK_LIBRARIES}")
    if(CP2K_LIBXS_INCLUDE_DIRS)
      if(CP2K_LIBXS_PREFIX)
        set_target_properties(
          cp2k::LIBXS::libxs
          PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                     "${CP2K_LIBXS_INCLUDE_DIRS};${CP2K_LIBXS_PREFIX}/include")
      else()
        set_target_properties(
          cp2k::LIBXS::libxs PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                        "${CP2K_LIBXS_INCLUDE_DIRS}")
      endif()
    endif()
  endif()
endif()

if(NOT TARGET cp2k::LIBXS)
  add_library(cp2k::LIBXS INTERFACE IMPORTED)
  target_link_libraries(cp2k::LIBXS INTERFACE cp2k::LIBXS::libxs)
endif()

mark_as_advanced(CP2K_LIBXS_INCLUDE_DIRS CP2K_LIBXS_LIBRARY_DIRS
                 CP2K_LIBXS_LINK_LIBRARIES)
