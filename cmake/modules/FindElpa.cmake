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

cp2k_set_default_paths(ELPA "Elpa")

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
  if(CP2K_ENABLE_ELPA_OPENMP_SUPPORT)
    pkg_search_module(CP2K_ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa_openmp)
  else()
    pkg_search_module(CP2K_ELPA REQUIRED IMPORTED_TARGET GLOBAL elpa)
  endif()
else()
  message(
    "ELPA : We need pkg-config to get all information about the elpa library")
endif()

find_package_handle_standard_args(Elpa "DEFAULT_MSG" CP2K_ELPA_LINK_LIBRARIES
                                  CP2K_ELPA_INCLUDE_DIRS)

if(CP2K_ELPA_FOUND)
  if(NOT TARGET cp2k::ELPA::elpa)
    add_library(cp2k::ELPA::elpa INTERFACE IMPORTED)
  endif()
  set_target_properties(
    cp2k::ELPA::elpa PROPERTIES INTERFACE_LINK_LIBRARIES
                                "${CP2K_ELPA_LINK_LIBRARIES}")
  if(CP2K_ELPA_INCLUDE_DIRS)
    set_target_properties(
      cp2k::ELPA::elpa
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                 "${CP2K_ELPA_INCLUDE_DIRS};${CP2K_ELPA_INCLUDE_DIRS}/modules")
  endif()
endif()

mark_as_advanced(CP2K_ELPA_FOUND CP2K_ELPA_LINK_LIBRARIES
                 CP2K_ELPA_INCLUDE_DIRS)
