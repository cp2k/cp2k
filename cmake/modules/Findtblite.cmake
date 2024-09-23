#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(TBLITE "tblite")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_TBLITE QUIET IMPORTED_TARGET GLOBAL tblite)
  pkg_check_modules(CP2K_DFTD3  QUIET IMPORTED_TARGET GLOBAL s-dftd3)
  pkg_check_modules(CP2K_TOML   QUIET IMPORTED_TARGET GLOBAL toml-f)
endif()

if(NOT CP2K_TBLITE_FOUND)
  cp2k_find_libraries(TBLITE tblite)
endif()
if(NOT CP2K_DFTD3_FOUND)
  cp2k_find_libraries(DFTD3 s-dftd3)
endif()
if(NOT CP2K_TOML_FOUND)
  cp2k_find_libraries(TOML toml-f)
endif()

if(NOT CP2K_TBLITE_INCLUDE_DIRS)
  cp2k_include_dirs(TBLITE "tblite.h")
endif()
if(NOT CP2K_DFTD3_INCLUDE_DIRS)
  cp2k_include_dirs(TBLITE "s-dftd3.h;dftd3.h")
endif()

if(CP2K_TBLITE_INCLUDE_DIRS)
  if(CP2K_DFTD3_INCLUDE_DIRS)
    find_package_handle_standard_args(tblite DEFAULT_MSG 
        CP2K_TBLITE_INCLUDE_DIRS CP2K_DFTD3_INCLUDE_DIRS
        CP2K_TBLITE_LINK_LIBRARIES CP2K_DFTD3_LINK_LIBRARIES
        CP2K_TOML_LINK_LIBRARIES)
  else()
    find_package_handle_standard_args(tblite DEFAULT_MSG 
        CP2K_TBLITE_INCLUDE_DIRS
        CP2K_TBLITE_LINK_LIBRARIES CP2K_DFTD3_LINK_LIBRARIES
        CP2K_TOML_LINK_LIBRARIES)
  endif()
else()
  if(CP2K_DFTD3_INCLUDE_DIRS)
    find_package_handle_standard_args(tblite DEFAULT_MSG 
        CP2K_DFTD3_INCLUDE_DIRS
        CP2K_TBLITE_LINK_LIBRARIES CP2K_DFTD3_LINK_LIBRARIES
        CP2K_TOML_LINK_LIBRARIES)
  else()
    find_package_handle_standard_args(tblite DEFAULT_MSG 
        CP2K_TBLITE_LINK_LIBRARIES CP2K_DFTD3_LINK_LIBRARIES
        CP2K_TOML_LINK_LIBRARIES)
  endif()
endif()

if(CP2K_TBLITE_FOUND)
  if(NOT TARGET cp2k::tblite::tblite)
    add_library(cp2k::tblite::tblite INTERFACE IMPORTED)
    if(NOT TARGET cp2k::tblite::s-dftd3)
      add_library(cp2k::tblite::s-dftd3 INTERFACE IMPORTED)
    endif()
    if(NOT TARGET cp2k::tblite::toml-f)
      add_library(cp2k::tblite::toml-f INTERFACE IMPORTED)
    endif()
  endif()
  
  set_target_properties(
     cp2k::tblite::tblite PROPERTIES INTERFACE_LINK_LIBRARIES "${CP2K_TBLITE_LINK_LIBRARIES}")
  set_target_properties(
     cp2k::tblite::s-dftd3 PROPERTIES INTERFACE_LINK_LIBRARIES "${CP2K_DFTD3_LINK_LIBRARIES}")
    set_target_properties(
     cp2k::tblite::toml-f PROPERTIES INTERFACE_LINK_LIBRARIES "${CP2K_TOML_LINK_LIBRARIES}")

    if(CP2K_TBLITE_INCLUDE_DIRS)
      set_target_properties(
        cp2k::tblite::tblite PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_TBLITE_INCLUDE_DIRS}" )
    endif()
    if(CP2K_TBLITE_INCLUDE_DIRS)
      set_target_properties(
        cp2k::tblite::s-dftd3 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_DFTD3_INCLUDE_DIRS}" )
    endif()             
endif()

mark_as_advanced(CP2K_TOML_LINK_LIBRARIES
                 CP2K_DFTD3_INCLUDE_DIRS  CP2K_DFTD3_LINK_LIBRARIES
                 CP2K_TBLITE_INCLUDE_DIRS CP2K_TBLITE_LINK_LIBRARIES
                 CP2K_TBLITE_FOUND)
