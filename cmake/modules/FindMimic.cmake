#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!





include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(PkgConfig REQUIRED)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_MIMIC IMPORTED_TARGET GLOBAL mimiccommf mimiccomm)
endif()

if(NOT CP2K_MIMIC_FOUND)
  cp2k_set_default_paths(MIMIC "Mimic")
  cp2k_find_libraries(MIMIC mimiccommf)
  cp2k_find_libraries(MIMICc mimiccomm)
endif()

if(CP2K_MIMIC_FOUND AND CP2K_MIMICc_FOUND)
  set(CP2K_MIMIC_LINK_LIBRARIES "${CP2K_MIMIC_LIBRARIES};${CP2K_MIMICc_LIBRARIES}")
endif()

if(NOT CP2K_MIMIC_INCLUDE_DIRS)
  cp2k_include_dirs(MIMIC "mcl.mod")
endif()

find_file(CP2K_MIMIC_MOD_FILE NAMES "mcl.mod" PATHS "${CP2K_MIMIC_ROOT}/include")

if(NOT CP2K_MIMIC_MOD_FILE)
  message(FATAL_ERROR "MIMIC : Fortran support in MCL is missing")
endif()

if(CP2K_MIMIC_INCLUDE_DIRS)
  find_package_handle_standard_args(Mimic DEFAULT_MSG CP2K_MIMIC_FOUND CP2K_MIMIC_LINK_LIBRARIES
                                                      CP2K_MIMIC_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(Mimic DEFAULT_MSG CP2K_MIMIC_FOUND CP2K_MIMIC_LINK_LIBRARIES)
endif()

if(CP2K_MIMIC_FOUND)
  if(NOT TARGET CP2K::MIMIC::mimiccommf)
    add_library(CP2K::MIMIC::mimiccommf INTERFACE IMPORTED)
    add_library(CP2K::MIMIC::mimiccomm INTERFACE IMPORTED)
  endif()

  if(CP2K_MIMIC_INCLUDE_DIRS)
    set_target_properties(
      CP2K::MIMIC::mimiccommf 
      CP2K::MIMIC::mimiccomm 
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_MIMIC_INCLUDE_DIRS}"
    )
  endif()
  target_link_libraries(
    CP2K::MIMIC::mimiccommf
    INTERFACE "${CP2K_MIMIC_LINK_LIBRARIES}"
  )
  target_link_libraries(
    CP2K::MIMIC::mimiccomm
    INTERFACE "${CP2K_MIMIC_LINK_LIBRARIES}"
  )
endif()

mark_as_advanced(CP2K_MIMIC_FOUND CP2K_MIMIC_LINK_LIBRARIES
                 CP2K_MIMIC_INCLUDE_DIRS)
