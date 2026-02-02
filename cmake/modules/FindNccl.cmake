#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(NCCL "Nccl")

cp2k_find_libraries(NCCL "nccl")
cp2k_include_dirs(NCCL "nccl.h")

find_package_handle_standard_args(Nccl DEFAULT_MSG CP2K_NCCL_LINK_LIBRARIES
                                  CP2K_NCCL_INCLUDE_DIRS)

if(CP2K_NCCL_FOUND AND NOT TARGET cp2k::NCCL::nccl)
  add_library(cp2k::NCCL::nccl INTERFACE IMPORTED)
  set_target_properties(
    cp2k::NCCL::nccl PROPERTIES INTERFACE_LINK_LIBRARIES
                                "${CP2K_NCCL_LINK_LIBRARIES}")
  set_target_properties(
    cp2k::NCCL::nccl PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                "${CP2K_NCCL_INCLUDE_DIRS}")
endif()

mark_as_advanced(CP2K_NCCL_LINK_LIBRARIES)
mark_as_advanced(CP2K_NCCL_INCLUDE_DIRS)
mark_as_advanced(CP2K_NCCL_FOUND)
