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

find_package(ptscotch)

cp2k_set_default_paths(PEXSI "PEXSI")

cp2k_find_libraries(PEXSI "pexsi")
cp2k_include_dirs(PEXSI "pexsi.hpp")

find_file(CP2K_PEXSI_MOD_FILE NAMES "f_ppexsi_interface.mod" PATCHS
                                    "${CP2K_PEXSI_PREFIX}/include")

if(NOT CP2K_PEXSI_MOD_FILE)
  message(
    FATAL_ERROR
      "The pexsi library needs to be compiled with fortran support. Either recompile pexsi or disable it."
  )
endif()

find_package_handle_standard_args(PEXSI DEFAULT_MSG CP2K_PEXSI_INCLUDE_DIRS
                                  CP2K_PEXSI_LINK_LIBRARIES)

if(CP2K_PEXSI_FOUND AND NOT TARGET cp2k::PEXSI::pexsi)
  add_library(cp2k::PEXSI::pexsi INTERFACE IMPORTED)
  set_target_properties(cp2k::PEXSI PROPERTIES INTERFACE_LINK_LIBRARIES
                                               "${CP2K_PEXSI_LINK_LIBRARIES}")
  if(DEFINED CP2K_PEXSI_INCLUDE_DIRS)
    set_target_properties(cp2k::PEXSI PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                 "${CP2K_PEXSI_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_PEXSI_LINK_LIBRARIES CP2K_PEXSI_INCLUDE_DIRS
                 CP2K_PEXSI_FOUND)
