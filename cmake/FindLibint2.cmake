include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(PkgConfig REQUIRED)
if(PKG_CONFIG_FOUND)
    pkg_check_modules(CP2K_LIBINT2 QUIET libint2)
endif()

if (NOT CP2K_LIBINT2_FOUND)
  cp2k_find_libraries(LIBINT2 int2)
  cp2k_include_dirs(LIBINT2 "libint2.h;libint2/atom.h")
endif()

find_package_handle_standard_args(
  Libint2
  CP2K_LIBINT2_FOUND
  CP2K_LIBINT2_INCLUDE_DIRS
  CP2K_LIBINT2_LINK_LIBRARIES)

if(CP2K_LIBINT2_FOUND AND NOT TARGET CP2K_Libint2::int2)
  add_library(CP2K_Libint2::int2 INTERFACE IMPORTED)
  set_target_properties(CP2K_Libint2::int2 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
    "${CP2K_LIBINT2_INCLUDE_DIRS}")
  set_target_properties(CP2K_Libint2::int2 PROPERTIES INTERFACE_LINK_LIBRARIES
    ${CP2K_LIBINT2_LINK_LIBRARIES})
endif()

mark_as_advanced(CP2K_LIBINT2_FOUND CP2K_LIBINT2_LINK_LIBRARIES CP2K_LIBINT2_INCLUDE_DIRS)
