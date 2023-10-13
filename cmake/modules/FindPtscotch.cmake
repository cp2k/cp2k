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

cp2k_set_default_paths(PTSCOTCH "Ptscotch")

find_package(Parmetis REQUIRED)
find_package(Threads REQUIRED)
find_package(MPI REQUIRED)

# look for libraries

foreach(
  _lib
  ptscotchparmetis
  ptscotch
  ptscotcherr
  scotchmetis
  scotch
  scotcherr
  ptesmumps)
  string(TOUPPER "${_lib}" _lib_up)
  cp2k_find_libraries("${_lib_up}" ${_lib})
endforeach()

# search for include files
cp2k_include_dirs(PTSCOTCH
                  "ptscotch.h openmpi/include/ptscotch.h ptsctoch/ptscotch.h")

# check that PTSCOTCH has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Ptscotch
  DEFAULT_MSG
  CP2K_PTSCOTCH_LINK_LIBRARIES
  CP2K_PTSCOTCHPARMETIS_LINK_LIBRARIES
  CP2K_PTSCOTCHERR_LINK_LIBRARIES
  CP2K_SCOTCHMETIS_LINK_LIBRARIES
  CP2K_SCOTCH_LINK_LIBRARIES
  CP2K_SCOTCHERR_LINK_LIBRARIES
  CP2K_PTESMUMPS_LINK_LIBRARIES)

if(CP2K_PTSCOTCH_FOUND AND NOT TARGET cp2k::ptscotch::ptscotch)
  add_library(cp2k::ptscotch::ptscotch INTERFACE IMPORTED)
  set_target_properties(
    cp2k::ptscotch::ptscotch
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${PTSCOTCH_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES
      "${CP2K_PTSCOTCH_LINK_LIBRARIES};
  ${CP2K_PTSCOTCHPARMETIS_LINK_LIBRARIES};
  ${CP2K_PTSCOTCHERR_LINK_LIBRARIES};
  ${CP2K_SCOTCHMETIS_LINK_LIBRARIES};
  ${CP2K_SCOTCH_LINK_LIBRARIES};
  ${CP2K_SCOTCHERR_LINK_LIBRARIES};
  ${CP2K_PTESMUMPS_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_PTSCOTCH_FOUND)
mark_as_advanced(CP2K_PTSCOTCH_LIBRARIES)
mark_as_advanced(CP2K_PTSCOTCH_INCLUDE_DIRS)
