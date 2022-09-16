# find spglib if in non-standard location set environment variabled `SPG_DIR` to
# the root directory

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(PkgConfig)

cp2k_set_default_paths(LIBSPG)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBSPG IMPORTED_TARGET spglib)
endif()

if(NOT CP2K_LIBSPG_FOUND)
  find_library(
    CP2K_LIBSPG_LIBRARIES
    NAMES "symspg"
    PATHS "${CP2K_LIBSPG_PREFIX}"
    HINTS "${CP2K_LIBSPG_PREFIX}")
  if (CP2K_LIBSPG_LIBRARIES)
    set(CP2K_LIBSPG_LINK_LIBRARIES ${CP2K_LIBSPG_LIBRARIES})
    set(CP2K_LIBSPG_FOUND TRUE)
  endif()
else()
  set(CP2K_LIBSPG_LINK_LIBRARIES ${CP2K_LIBSPG_LINK_LIBRARIES})
endif()

if (NOT DEFINED CP2K_LIBSPG_INCLUDE_DIRS)
  find_path(CP2K_LIBSPG_INCLUDE_DIRS
    NAMES "spglib.h"
    PATHS "${CP2K_LIBSPG_PREFIX}"
    HINTS "${CP2K_LIBSPG_PREFIX}"
    PATH_SUFFIXES "spglib" "include/spglib" 
    NO_DEFAULT_PATH)
endif()

if(CP2K_LIBSPG_INCLUDE_DIRS)
  find_package_handle_standard_args(LibSPG DEFAULT_MSG CP2K_LIBSPG_LINK_LIBRARIES
    CP2K_LIBSPG_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibSPG DEFAULT_MSG CP2K_LIBSPG_LINK_LIBRARIES)
endif()

if(LIBSPG_FOUND AND NOT TARGET CP2K_LIBSPG::libspg)
  add_library(CP2K_LIBSPG::libspg INTERFACE IMPORTED)
  set_target_properties(CP2K_LIBSPG::libspg PROPERTIES INTERFACE_LINK_LIBRARIES
    "${CP2K_LIBSPG_LINK_LIBRARIES}")
  if(CP2K_LIBSPG_INCLUDE_DIRS)
    set_target_properties(
      CP2K_LIBSPG::libspg PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
      "${CP2K_LIBSPG_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBSPG_LINK_LIBRARIES)
mark_as_advanced(CP2K_LIBSPG_INCLUDE_DIRS)
mark_as_advanced(CP2K_LIBSPG_FOUND)
