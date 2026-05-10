#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(libfci CONFIG QUIET)

if(libfci_FOUND AND TARGET libfci::fci)
  set(LibFCI_FOUND TRUE)
  set(CP2K_LIBFCI_FOUND TRUE)
  set(CP2K_LIBFCI_LINK_LIBRARIES libfci::fci)
  get_target_property(CP2K_LIBFCI_INCLUDE_DIRS libfci::fci
                      INTERFACE_INCLUDE_DIRECTORIES)
else()
  find_package(PkgConfig QUIET)
  cp2k_set_default_paths(LIBFCI "libfci")

  if(PKG_CONFIG_FOUND)
    pkg_check_modules(CP2K_LIBFCI IMPORTED_TARGET GLOBAL libfci)
  endif()

  if(NOT CP2K_LIBFCI_FOUND)
    cp2k_find_libraries(LIBFCI "fci")
    cp2k_include_dirs(LIBFCI "libfci.h")
  endif()

  find_package_handle_standard_args(LibFCI DEFAULT_MSG CP2K_LIBFCI_INCLUDE_DIRS
                                    CP2K_LIBFCI_LINK_LIBRARIES)
endif()

if(CP2K_LIBFCI_FOUND)
  if(NOT TARGET cp2k::libfci::fci)
    add_library(cp2k::libfci::fci INTERFACE IMPORTED)
  endif()

  if(TARGET libfci::fci)
    target_link_libraries(cp2k::libfci::fci INTERFACE libfci::fci)
  elseif(TARGET PkgConfig::CP2K_LIBFCI)
    target_link_libraries(cp2k::libfci::fci INTERFACE PkgConfig::CP2K_LIBFCI)
  else()
    set_target_properties(
      cp2k::libfci::fci PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                   "${CP2K_LIBFCI_INCLUDE_DIRS}")
    target_link_libraries(cp2k::libfci::fci
                          INTERFACE ${CP2K_LIBFCI_LINK_LIBRARIES})
  endif()
endif()

mark_as_advanced(CP2K_LIBFCI_FOUND CP2K_LIBFCI_INCLUDE_DIRS
                 CP2K_LIBFCI_LINK_LIBRARIES)
