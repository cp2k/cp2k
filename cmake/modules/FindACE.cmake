#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(ACE "ace")

cp2k_include_dirs(
  ACE "ace/ace_couplings.h;ace-evaluator/ace_types.h;yaml-cpp/yaml.h")
find_package_handle_standard_args(ACE DEFAULT_MSG CP2K_ACE_INCLUDE_DIRS)

cp2k_find_libraries(ACE "pace")
find_package_handle_standard_args(ACE DEFAULT_MSG CP2K_ACE_LINK_LIBRARIES)

cp2k_find_libraries(ACE_YAML "yaml-cpp-pace")
find_package_handle_standard_args(ACE DEFAULT_MSG CP2K_ACE_YAML_LINK_LIBRARIES)

cp2k_find_libraries(ACE_CNPY "cnpy")
find_package_handle_standard_args(ACE DEFAULT_MSG CP2K_ACE_CNPY_LINK_LIBRARIES)

if(CP2K_ACE_FOUND)
  if(NOT TARGET ACE::pace)
    add_library(ACE::pace INTERFACE IMPORTED)
  endif()
  set_target_properties(
    ACE::pace
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_ACE_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_ACE_LINK_LIBRARIES}")
endif()

if(CP2K_ACE_YAML_FOUND)
  if(NOT TARGET ACE::yaml-cpp-pace)
    add_library(ACE::yaml-cpp-pace INTERFACE IMPORTED)
  endif()
  set_target_properties(
    ACE::yaml-cpp-pace PROPERTIES INTERFACE_LINK_LIBRARIES
                                  "${CP2K_ACE_YAML_LINK_LIBRARIES}")
endif()

if(CP2K_ACE_CNPY_FOUND)
  if(NOT TARGET ACE::cnpy)
    add_library(ACE::cnpy INTERFACE IMPORTED)
  endif()
  set_target_properties(ACE::cnpy PROPERTIES INTERFACE_LINK_LIBRARIES
                                             "${CP2K_ACE_CNPY_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_ACE_FOUND CP2K_ACE_INCLUDE_DIRS CP2K_ACE_LINK_LIBRARIES
                 CP2K_ACE_YAML_LINK_LIBRARIES CP2K_ACE_CNPY_LINK_LIBRARIES)
