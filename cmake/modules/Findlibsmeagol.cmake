#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2025- ETH Zurich
#
# author: Rocco Meli

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_find_libraries(LIBSMEAGOL "smeagol")
cp2k_include_dirs(LIBSMEAGOL "negfmod.mod")

if(CP2K_LIBSMEAGOL_INCLUDE_DIRS)
  find_package_handle_standard_args(
    libsmeagol DEFAULT_MSG CP2K_LIBSMEAGOL_LINK_LIBRARIES
    CP2K_LIBSMEAGOL_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(libsmeagol DEFAULT_MSG
                                    CP2K_LIBSMEAGOL_LINK_LIBRARIES)
endif()

if(NOT TARGET cp2k::libsmeagol::smeagol)
  add_library(cp2k::libsmeagol::smeagol INTERFACE IMPORTED)
  set_target_properties(
    cp2k::libsmeagol::smeagol PROPERTIES INTERFACE_LINK_LIBRARIES
                                         "${CP2K_LIBSMEAGOL_LINK_LIBRARIES}")
  if(CP2K_LIBSMEAGOL_INCLUDE_DIRS)
    set_target_properties(
      cp2k::libsmeagol::smeagol PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                           "${CP2K_LIBSMEAGOL_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBSMEAGOL_ROOT CP2K_LIBSMEAGOL_INCLUDE_DIRS
                 CP2K_LIBSMEAGOL_LINK_LIBRARIES CP2K_LIBSMEAGOL_LIBRARIES)
