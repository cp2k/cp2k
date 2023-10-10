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

cp2k_set_default_paths(METIS "Metis")

cp2k_find_libraries(FLEXIBLAS "metis")
cp2k_include_dirs(FFTW3 "metis.h")

# check that METIS has been found
# ---------------------------------
find_package_handle_standard_args(Metis DEFAULT_MSG CP2K_METIS_LINK_LIBRARIES
                                  CP2K_METIS_INCLUDE_DIRS CP2K_METIS_FOUND)

if(CP2K_METIS_FOUND AND NOT TARGET cp2k::metis::metis)
  add_library(cp2k::metis::metis INTERFACE IMPORTED)
  set_target_properties(
    cp2k::metis::metis
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_METIS_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_METIS_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_METIS_LINK_LIBRARIES)
mark_as_advanced(CP2K_METIS_INCLUDE_DIRS)
mark_as_advanced(CP2K_METIS_FOUND)
