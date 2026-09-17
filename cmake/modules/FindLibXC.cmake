#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(LIBXC "libxc")

if(PKG_CONFIG_FOUND)
  # Use the static resolution so the C library is pulled in explicitly
  # (Requires.private libxc), then turn the archive names into full paths for
  # linking.
  pkg_check_modules(CP2K_LIBXC IMPORTED_TARGET GLOBAL libxcf03)
  if(CP2K_LIBXC_FOUND)
    unset(CP2K_LIBXC_LINK_LIBRARIES)
    foreach(_lib IN LISTS CP2K_LIBXC_STATIC_LIBRARIES)
      unset(_libxc_library CACHE)
      find_library(
        _libxc_library
        NAMES ${_lib}
        PATHS ${CP2K_LIBXC_STATIC_LIBRARY_DIRS})
      if(_libxc_library)
        list(APPEND CP2K_LIBXC_LINK_LIBRARIES ${_libxc_library})
      else()
        list(APPEND CP2K_LIBXC_LINK_LIBRARIES ${_lib})
      endif()
    endforeach()
    set(CP2K_LIBXC_INCLUDE_DIRS "${CP2K_LIBXC_INCLUDEDIR}")
  endif()
endif()

if(NOT CP2K_LIBXC_FOUND)
  cp2k_find_libraries(LIBXC "xcf03")
endif()

if(NOT CP2K_LIBXC_INCLUDE_DIRS)
  cp2k_include_dirs(LIBXC "xc.h")
endif()

find_package_handle_standard_args(LibXC DEFAULT_MSG CP2K_LIBXC_INCLUDE_DIRS
                                  CP2K_LIBXC_LINK_LIBRARIES)

if(CP2K_LIBXC_FOUND)
  if(NOT TARGET Libxc::xcf03)
    add_library(Libxc::xcf03 INTERFACE IMPORTED)
  endif()
  set_target_properties(
    Libxc::xcf03
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXC_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_LIBXC_LINK_LIBRARIES}")

  if(NOT TARGET Libxc::xc)
    add_library(Libxc::xc INTERFACE IMPORTED)
  endif()
  set_target_properties(
    Libxc::xc
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBXC_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_LIBXC_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_LIBXC_FOUND CP2K_LIBXC_INCLUDE_DIRS
                 CP2K_LIBXC_LINK_LIBRARIES)
