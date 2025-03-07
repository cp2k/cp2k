#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(CAL "Cal")

cp2k_find_libraries(CAL "cal")
cp2k_include_dirs(CAL "cal.h")

find_package_handle_standard_args(Cal DEFAULT_MSG CP2K_CAL_LINK_LIBRARIES
                                  CP2K_CAL_INCLUDE_DIRS)

if(CP2K_CAL_FOUND AND NOT TARGET cp2k::CAL::cal)
  add_library(cp2k::CAL::cal INTERFACE IMPORTED)
  set_target_properties(cp2k::CAL::cal PROPERTIES INTERFACE_LINK_LIBRARIES
                                                  "${CP2K_CAL_LINK_LIBRARIES}")
  set_target_properties(cp2k::CAL::cal PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                  "${CP2K_CAL_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Cal required by CuSolverMP")
endif()

mark_as_advanced(CP2K_CAL_LINK_LIBRARIES)
mark_as_advanced(CP2K_CAL_INCLUDE_DIRS)
mark_as_advanced(CP2K_CAL_FOUND)
