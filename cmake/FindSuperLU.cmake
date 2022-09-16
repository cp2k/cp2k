include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(PkgConfig)
cp2k_set_default_paths(SUPERLU)
pkg_search_module(CP2K_SUPERLU "superlu superlu_dist")

if (NOT CP2K_SUPERLU_FOUND)
  cp2k_find_libraries(SUPERLU "superlu")
  cp2k_include_dirs(SUPERLU "supermatrix.h;SuperLU/supermatrix.h;superlu/supermatrix.h")
endif()

find_package_handle_standard_args(SuperLU DEFAULT_MSG CP2K_SUPERLU_INCLUDE_DIRS
  CP2K_SUPERLU_LINK_LIBRARIES)

if(CP2K_SUPERLU_FOUND AND NOT TARGET CP2K_superlu::superlu)
  add_library(CP2K_superlu::superlu INTERFACE IMPORTED)
  set_target_properties(
    CP2K_superlu::superlu
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_SUPERLU_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${CP2K_SUPERLU_LINK_LIBRARIES}")
endif()
