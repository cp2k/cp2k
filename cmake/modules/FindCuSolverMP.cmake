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

find_package(Cal Required)
cp2k_set_default_paths(CUSOLVER_MP "CUSOLVER_MP")
cp2k_find_libraries(CUSOLVER_MP "cusolverMp")
cp2k_include_dirs(CUSOLVER_MP "cusolverMp.h")

find_package_handle_standard_args(
  CuSolverMP DEFAULT_MSG CP2K_CUSOLVER_MP_LINK_LIBRARIES
  CP2K_CUSOLVER_MP_INCLUDE_DIRS)

if(CP2K_CUSOLVER_MP_FOUND AND NOT TARGET cp2k::CUSOLVER_MP::cusolver_mp)
  add_library(cp2k::CUSOLVER_MP::cusolver_mp INTERFACE IMPORTED)
  set_target_properties(
    cp2k::CUSOLVER_MP::cusolver_mp
    PROPERTIES INTERFACE_LINK_LIBRARIES
               "${CP2K_CUSOLVER_MP_LINK_LIBRARIES};cp2k::CAL::cal")
  set_target_properties(
    cp2k::CUSOLVER_MP::cusolver_mp
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_CUSOLVER_MP_INCLUDE_DIRS}")
endif()

mark_as_advanced(CP2K_CUSOLVER_MP_LINK_LIBRARIES)
mark_as_advanced(CP2K_CUSOLVER_MP_INCLUDE_DIRS)
mark_as_advanced(CP2K_CUSOLVER_MP_FOUND)
