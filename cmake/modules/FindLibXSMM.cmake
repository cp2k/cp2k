#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)
find_package(PkgConfig REQUIRED)

cp2k_set_default_paths(LIBXSMM "LibXSMM")
set(CP2K_LIBXSMMEXT_ROOT "${CP2K_LIBXSMM_PREFIX}")
set(CP2K_LIBXSMMF_ROOT "${CP2K_LIBXSMM_PREFIX}")
set(CP2K_LIBXSMMNOBLAS_ROOT "${CP2K_LIBXSMM_PREFIX}")

if(PKG_CONFIG_FOUND)
  foreach(__lib libxsmm libxsmmf libxsmmext libxsmmnoblas)
    string(TOUPPER "${__lib}" __lib_search_up)
    pkg_check_modules(CP2K_${__lib_search_up} IMPORTED_TARGET GLOBAL ${__lib})
    # need to do it twice because of dbcsr build option
    pkg_check_modules(${__lib_search_up} QUIET IMPORTED_TARGET GLOBAL ${__lib})
  endforeach()
endif()

if(NOT CP2K_LIBXSMM_FOUND)
  # Reset after pkg_check_modules side effects
  foreach(__lib xsmm xsmmf xsmmext xsmmnoblas)
    string(TOUPPER "LIB${__lib}" __lib_search_up)
    if(NOT CP2K_${__lib_search_up}_FOUND)
      cp2k_find_libraries(${__lib_search_up} ${__lib})
    endif()
  endforeach()
endif()

if(NOT CP2K_LIBXSMM_INCLUDE_DIRS)
  cp2k_include_dirs(LIBXSMM "libxsmm.h;include/libxsmm.h")
endif()

if(CP2K_LIBXSMM_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibXSMM
    DEFAULT_MSG
    CP2K_LIBXSMM_INCLUDE_DIRS
    CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
    CP2K_LIBXSMMEXT_LINK_LIBRARIES
    CP2K_LIBXSMMF_LINK_LIBRARIES
    CP2K_LIBXSMM_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(
    LibXSMM DEFAULT_MSG CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
    CP2K_LIBXSMMEXT_LINK_LIBRARIES CP2K_LIBXSMMF_LINK_LIBRARIES
    CP2K_LIBXSMM_LINK_LIBRARIES)
endif()

if(NOT TARGET cp2k::LibXSMM::libxsmm)
  foreach(__lib libxsmm libxsmmf libxsmmext libxsmmnoblas)
    string(TOUPPER "CP2K_${__lib}" __lib_search_up)

    if(${__lib_search_up}_FOUND AND NOT TARGET cp2k::LibXSMM::${__lib})
      add_library(cp2k::LibXSMM::${__lib} INTERFACE IMPORTED)
      target_link_directories(cp2k::LibXSMM::${__lib} INTERFACE
                              ${${__lib_search_up}_LIBRARY_DIRS})
    endif()

    set_target_properties(
      cp2k::LibXSMM::${__lib} PROPERTIES INTERFACE_LINK_LIBRARIES
                                         "${${__lib_search_up}_LINK_LIBRARIES}")

    if(CP2K_LIBXSMM_INCLUDE_DIRS)
      set_target_properties(
        cp2k::LibXSMM::${__lib}
        PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                   "${CP2K_LIBXSMM_INCLUDE_DIRS};${CP2K_LIBXSMM_PREFIX}/include"
      )
    endif()
  endforeach()
endif()

mark_as_advanced(
  CP2K_LIBXSMM_INCLUDE_DIRS CP2K_LIBXSMM_LIBRARY_DIRS
  CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES CP2K_LIBXSMMEXT_LINK_LIBRARIES
  CP2K_LIBXSMMF_LINK_LIBRARIES CP2K_LIBXSMM_LINK_LIBRARIES)
