#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(ACE "ace")

cp2k_include_dirs(
  ACE "ace/ace_couplings.h;ace-evaluator/ace_types.h;yaml-cpp/yaml.h")
cp2k_find_libraries(ACE "pace")
cp2k_find_libraries(ACE_YAML "yaml-cpp-pace")
cp2k_find_libraries(ACE_CNPY "cnpy")

find_package_handle_standard_args(
  ACE DEFAULT_MSG CP2K_ACE_LINK_LIBRARIES CP2K_ACE_CNPY_LINK_LIBRARIES
  CP2K_ACE_INCLUDE_DIRS CP2K_ACE_YAML_LINK_LIBRARIES)

if(CP2K_ACE_FOUND)
  if(NOT TARGET cp2k::ACE)
    add_library(cp2k::ACE INTERFACE IMPORTED)
  endif()
  set_target_properties(
    cp2k::ACE
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CP2K_ACE_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES
      "${CP2K_ACE_LINK_LIBRARIES};${CP2K_ACE_YAML_LINK_LIBRARIES};${CP2K_ACE_CNPY_LINK_LIBRARIES}"
  )
endif()

mark_as_advanced(CP2K_ACE_FOUND CP2K_ACE_INCLUDE_DIRS CP2K_ACE_LINK_LIBRARIES
                 CP2K_ACE_YAML_LINK_LIBRARIES CP2K_ACE_CNPY_LINK_LIBRARIES)
