#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(TREXIO "trexio")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_TREXIO IMPORTED_TARGET GLOBAL trexio)
endif()

if(NOT CP2K_TREXIO_FOUND)
  cp2k_find_libraries(TREXIO "trexio")
  cp2k_include_dirs(TREXIO "trexio.h")
endif()

find_package_handle_standard_args(TrexIO DEFAULT_MSG CP2K_TREXIO_INCLUDE_DIRS
                                  CP2K_TREXIO_LINK_LIBRARIES)

if(CP2K_TREXIO_FOUND)
  if(NOT TARGET cp2k::trexio::trexio)
    add_library(cp2k::trexio::trexio INTERFACE IMPORTED)
  endif()
  set_target_properties(
    cp2k::trexio::trexio PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                    "${CP2K_TREXIO_INCLUDE_DIRS}")
  target_link_libraries(cp2k::trexio::trexio
                        INTERFACE ${CP2K_TREXIO_LINK_LIBRARIES})
endif()

mark_as_advanced(CP2K_TREXIO_FOUND CP2K_TREXIO_INCLUDE_DIRS
                 CP2K_TREXIO_LINK_LIBRARIES)
