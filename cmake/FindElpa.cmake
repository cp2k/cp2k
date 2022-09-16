# find Elpa via pkg-config or easybuild

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(ELPA)

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
  pkg_search_module(
    CP2K_ELPA
    elpa
    elpa_openmp
    elpa-openmp-2021.11.002
    elpa-2021.11.002)
else()
  message("ELPA : We need pkg-config to get all information about the elpa library")
endif()

if(NOT CP2K_ELPA_INCLUDE_DIRS)
  cp2k_include_dirs(ELPA "elpa.h")
endif()

find_package_handle_standard_args(Elpa "DEFAULT_MSG" CP2K_ELPA_LINK_LIBRARIES
                                  CP2K_ELPA_INCLUDE_DIRS)

if(CP2K_ELPA_FOUND AND NOT TARGET CP2K_elpa::elpa)
  add_library(CP2K_ELPA::elpa INTERFACE IMPORTED)
  set_target_properties(CP2K_ELPA::elpa PROPERTIES INTERFACE_LINK_LIBRARIES
    "${CP2K_ELPA_LINK_LIBRARIES}")
  if(CP2K_ELPA_INCLUDE_DIRS)
    set_target_properties(CP2K_ELPA::elpa PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
      "${CP2K_ELPA_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_ELPA_FOUND CP2K_ELPA_LINK_LIBRARIES CP2K_ELPA_INCLUDE_DIRS)
