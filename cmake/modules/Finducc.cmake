#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# authors : JVP find ucc needed CuSolverMP

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(UCC "ucc")

cp2k_find_libraries(UCC "ucc")
cp2k_find_libraries(UCX "ucs")

if(NOT CP2K_UCC_INCLUDE_DIRS)
  cp2k_include_dirs(UCC "ucc.h;ucc/api/ucc.h")
endif()

find_package_handle_standard_args(
  ucc DEFAULT_MSG CP2K_UCC_INCLUDE_DIRS CP2K_UCC_LINK_LIBRARIES
  CP2K_UCX_LINK_LIBRARIES)

if(CP2K_UCC_FOUND AND NOT TARGET cp2k::UCC::ucc)
  add_library(cp2k::UCC::ucc INTERFACE IMPORTED)
  if(CP2K_UCX_FOUND)
    set_target_properties(
      cp2k::UCC::ucc PROPERTIES INTERFACE_LINK_LIBRARIES
                                "${CP2K_UCC_LINK_LIBRARIES}")
  else()
    set_target_properties(
      cp2k::UCC::ucc
      PROPERTIES INTERFACE_LINK_LIBRARIES
                 "${CP2K_UCC_LINK_LIBRARIES};${CP2K_UCX_LINK_LIBRARIES}")
  endif()
  set_target_properties(cp2k::UCC::ucc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                  "${CP2K_UCC_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "ucc required by CuSolverMP")
endif()

mark_as_advanced(CP2K_UCX_LINK_LIBRARIES)
mark_as_advanced(CP2K_UCC_LINK_LIBRARIES)
mark_as_advanced(CP2K_UCC_INCLUDE_DIRS)
mark_as_advanced(CP2K_UCC_FOUND)
