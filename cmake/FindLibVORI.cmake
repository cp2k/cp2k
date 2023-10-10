#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(LIBVORI "LibVORI")
cp2k_find_libraries(LIBVORI vori)
# cp2k_include_dirs(LIBVORI )

if(CP2K_LIBVORI_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibVORI DEFAULT_MSG CP2K_LIBVORI_LINK_LIBRARIES CP2K_LIBVORI_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibVORI DEFAULT_MSG
                                    CP2K_LIBVORI_LINK_LIBRARIES)
endif()

if(NOT TARGET cp2k::VORI::vori)
  add_library(cp2k::VORI::vori INTERFACE IMPORTED)
  set_target_properties(
    cp2k::VORI::vori PROPERTIES INTERFACE_LINK_LIBRARIES
                                "${CP2K_LIBVORI_LINK_LIBRARIES}")
  if(CP2K_LIBVORI_INCLUDE_DIRS)
    set_target_properties(
      cp2k::VORI::vori PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                  "${CP2K_LIBVORI_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBVORI_ROOT CP2K_LIBVORI_INCLUDE_DIRS
                 CP2K_LIBVORI_LINK_LIBRARIES CP2K_LIBVORI_LIBRARIES)
