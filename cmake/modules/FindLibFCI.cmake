#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(LIBFCI "LibFCI")
cp2k_find_libraries(LIBFCI "fci")

if(CP2K_LIBFCI_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibFCI DEFAULT_MSG CP2K_LIBFCI_LINK_LIBRARIES CP2K_LIBFCI_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibFCI DEFAULT_MSG
                                    CP2K_LIBFCI_LINK_LIBRARIES)
endif()

if(LibFCI_FOUND)
  if(NOT TARGET libfci::fci)
    add_library(libfci::fci INTERFACE IMPORTED)
    set_target_properties(
      libfci::fci PROPERTIES INTERFACE_LINK_LIBRARIES
                             "${CP2K_LIBFCI_LINK_LIBRARIES}")
    if(CP2K_LIBFCI_INCLUDE_DIRS)
      set_target_properties(
        libfci::fci PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                               "${CP2K_LIBFCI_INCLUDE_DIRS}")
    endif()
  endif()
endif()

mark_as_advanced(CP2K_LIBFCI_ROOT CP2K_LIBFCI_INCLUDE_DIRS
                 CP2K_LIBFCI_LINK_LIBRARIES CP2K_LIBFCI_LIBRARIES)
