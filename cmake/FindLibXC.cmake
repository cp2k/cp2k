#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)
include(cp2k_utils)

cp2k_set_default_paths(LIBXC "LibXC")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBXC IMPORTED_TARGET GLOBAL libxcf90 libxcf03
                    libxc>=${LibXC_FIND_VERSION})
endif()

if(NOT CP2K_LIBXC_FOUND)
  # Revert pkg_check_modules side effects
  cp2k_set_default_paths(LIBXC "LibXC")
  foreach(_var xc xcf03 xcf90)
    string(TOUPPER LIB${_var} _var_up)
    cp2k_find_libraries(${_var_up} ${_var})
  endforeach()
endif()

if(CP2K_LIBXC_FOUND
   AND CP2K_LIBXCF90_FOUND
   AND CP2K_LIBXCF03_FOUND)
  set(CP2K_LIBXC_LINK_LIBRARIES
      "${CP2K_LIBXCF03_LIBRARIES};${CP2K_LIBXCF90_LIBRARIES};${CP2K_LIBXC_LIBRARIES}"
  )
endif()

if(NOT CP2K_LIBXC_INCLUDE_DIRS)
  cp2k_include_dirs(LIBXC "xc.h;libxc/xc.h")
endif()

if(CP2K_LIBXC_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibXC DEFAULT_MSG CP2K_LIBXC_FOUND CP2K_LIBXC_LINK_LIBRARIES
    CP2K_LIBXC_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibXC DEFAULT_MSG CP2K_LIBXC_FOUND
                                    CP2K_LIBXC_LINK_LIBRARIES)
endif()
if(CP2K_LIBXC_FOUND)
  if(NOT TARGET cp2k::Libxc::xc)
    add_library(cp2k::Libxc::xc INTERFACE IMPORTED)
  endif()

  if(CP2K_LIBXC_INCLUDE_DIRS)
    set_target_properties(
      cp2k::Libxc::xc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                 "${CP2K_LIBXC_INCLUDE_DIRS}")
  endif()
  target_link_libraries(cp2k::Libxc::xc INTERFACE ${CP2K_LIBXC_LINK_LIBRARIES})
endif()

mark_as_advanced(CP2K_LIBXC_FOUND CP2K_LIBXC_LINK_LIBRARIES
                 CP2K_LIBXC_INCLUDE_DIRS)
