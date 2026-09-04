#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Find libwignernj, which provides the exact evaluation of the Wigner symbols,
# Clebsch-Gordan coefficients and Gaunt coefficients.
#
# CP2K binds to the C interface directly through ISO_C_BINDING, so only the C
# library is needed; the Fortran interface library of libwignernj is not used.

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(PkgConfig)

cp2k_set_default_paths(LIBWIGNERNJ "Libwignernj")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBWIGNERNJ IMPORTED_TARGET GLOBAL libwignernj)
endif()

if(NOT CP2K_LIBWIGNERNJ_FOUND)
  cp2k_find_libraries(LIBWIGNERNJ wignernj)
endif()

if(NOT DEFINED CP2K_LIBWIGNERNJ_INCLUDE_DIRS)
  cp2k_include_dirs(LIBWIGNERNJ "wignernj.h")
endif()

find_package_handle_standard_args(Libwignernj DEFAULT_MSG
                                  CP2K_LIBWIGNERNJ_LINK_LIBRARIES)

if(CP2K_LIBWIGNERNJ_FOUND AND NOT TARGET cp2k::libwignernj::wignernj)
  add_library(cp2k::libwignernj::wignernj INTERFACE IMPORTED)
  set_target_properties(
    cp2k::libwignernj::wignernj PROPERTIES INTERFACE_LINK_LIBRARIES
                                           "${CP2K_LIBWIGNERNJ_LINK_LIBRARIES}")
  if(CP2K_LIBWIGNERNJ_INCLUDE_DIRS)
    set_target_properties(
      cp2k::libwignernj::wignernj PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                             "${CP2K_LIBWIGNERNJ_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBWIGNERNJ_ROOT CP2K_LIBWIGNERNJ_INCLUDE_DIRS
                 CP2K_LIBWIGNERNJ_LINK_LIBRARIES CP2K_LIBWIGNERNJ_FOUND)
