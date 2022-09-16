include(FindPackageHandleStandardArgs)
include(cp2k_utils)
find_package(PkgConfig)

cp2k_set_default_paths(PLUMED)

# First try with pkg
if(PKG_CONFIG_FOUND)
  # plumed has a pkg-config module
  pkg_check_module(CP2K_PLUMED "plumed")
endif()

if (NOT ${CP2K_PLUMED_FOUND})
  cp2k_find_libraries(PLUMED "plumed")
endif()

cp2k_include_dirs(PLUMED "plumed.h plumed/plumed.h")

if (CP2K_PLUMED_INCLUDE_DIRS)
find_package_handle_standard_args(Plumed DEFAULT_MSG CP2K_PLUMED_INCLUDE_DIRS
  CP2K_PLUMED_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Plumed DEFAULT_MSG 
    CP2K_PLUMED_LINK_LIBRARIES)
endif()

if(CP2K_PLUMED_FOUND AND NOT TARGET CP2K_plumed::plumed)
  add_library(CP2K_plumed::plumed INTERFACE IMPORTED)
  set_target_properties(
    CP2K_plumed::plumed
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_PLUMED_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${CP2K_PLUMED_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_PLUMED_LINK_LIBRARIES CP2K_PLUMED_INCLUDE_DIRS CP2K_PLUMED_FOUND)
