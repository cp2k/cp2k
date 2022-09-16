include(FindPackageHandleStandardArgs)
find_package(PkgConfig)
cp2k_set_default_paths(LIBQUIP)

pkg_search_module(CP2K_LIBQUIP quip)

if (CP2K_LIBQUIP_FOUND)
  cp2k_find_libraries(LIBQUIP "quip")
  cp2k_include_dirs(LIBQUIP "quip.h")
endif()
find_package_handle_standard_args(Quip DEFAULT_MSG CP2K_LIBQUIP_INCLUDE_DIRS
  CP2K_QUIP_LINK_LIBRARIES)

if(CP2K_LIBQUIP_FOUND AND NOT TARGET CP2K_quip::quip)
  add_library(CP2K_quip::quip INTERFACE IMPORTED)
  set_target_properties(
    CP2K_quip::quip PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBQUIP_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${CP2K_LIBQUIP_LINK_LIBRARIES}")
endif()
