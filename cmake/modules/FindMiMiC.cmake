#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

# ! Use PkgConfig to look for modules named mclf and mcl, LOOK for .pc files
find_package(PkgConfig REQUIRED)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_MIMIC IMPORTED_TARGET GLOBAL mclf mcl)
endif()

# ! Using CP2K_UTILS
if(NOT CP2K_MIMIC_FOUND)
  cp2k_set_default_paths(MIMIC "MiMiC")
  cp2k_find_libraries(MIMIC mclf)
  cp2k_find_libraries(MIMICc mcl)
endif()

if(NOT CP2K_MIMIC_FOUND)
  find_library(CP2K_MIMIC_LIBRARIES mclf PATH_SUFFIXES MiMiC)
  find_library(CP2K_MIMICc_LIBRARIES mcl PATH_SUFFIXES MiMiC)
  set(CP2K_MIMIC_FOUND True)
endif()

if(CP2K_MIMIC_FOUND)
  set(CP2K_MIMIC_LINK_LIBRARIES
      "${CP2K_MIMIC_LIBRARIES};${CP2K_MIMICc_LIBRARIES}")
endif()

if(NOT CP2K_MIMIC_INCLUDE_DIRS)
  cp2k_include_dirs(MIMIC "mcl.mod")
endif()

if(NOT CP2K_MIMIC_INCLUDE_DIRS)
  find_path(CP2K_MIMIC_INCLUDE_DIRS "mcl.mod" PATH_SUFFIXES MiMiC)
endif()

if(CP2K_MIMIC_INCLUDE_DIRS)
  find_package_handle_standard_args(
    MiMiC DEFAULT_MSG CP2K_MIMIC_FOUND CP2K_MIMIC_LINK_LIBRARIES
    CP2K_MIMIC_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(MiMiC DEFAULT_MSG CP2K_MIMIC_FOUND
                                    CP2K_MIMIC_LINK_LIBRARIES)
endif()

if(CP2K_MIMIC_FOUND)
  if(NOT TARGET CP2K::MIMIC::mclf)
    add_library(CP2K::MIMIC::mclf INTERFACE IMPORTED)
    add_library(CP2K::MIMIC::mcl INTERFACE IMPORTED)
  endif()

  if(CP2K_MIMIC_INCLUDE_DIRS)
    set_target_properties(
      CP2K::MIMIC::mclf CP2K::MIMIC::mcl
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_MIMIC_INCLUDE_DIRS}")
  endif()
  target_link_libraries(CP2K::MIMIC::mclf
                        INTERFACE "${CP2K_MIMIC_LINK_LIBRARIES}")
  target_link_libraries(CP2K::MIMIC::mcl
                        INTERFACE "${CP2K_MIMIC_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_MIMIC_FOUND CP2K_MIMIC_LINK_LIBRARIES
                 CP2K_MIMIC_INCLUDE_DIRS)
