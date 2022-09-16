
include(FindPackageHandleStandardArgs)
include(cp2k_utils)
include(FindPackageHandleStandardArgs)

find_package(ptscotch)

cp2k_set_default_paths(PEXSI)

cp2k_find_libraries(PEXSI pexsi)
cp2k_include_dirs(PEXSI "pexsi.h pexsi/pexsi.h")

find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_PEXSI_INCLUDE_DIRS
  CP2K_PEXSI_LINK_LIBRARIES)

if (CP2K_PEXSI_FOUND AND NOT TARGET CP2K_PEXSI::pexsi)
  add_library(CP2K_PEXSI::pexsi INTERFACE IMPORTED)
  set_target_properties(CP2K_PEXSI PROPERTIES INTERFACE_LINK_LIBRARIES
    "${CP2K_PEXSI_LINK_LIBRARIES}")
  if (DEFINED CP2K_PEXSI_INCLUDE_DIRS)
    set_target_properties(CP2K_PEXSI PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
      "${CP2K_PEXSI_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_PEXSI_LINK_LIBRARIES CP2K_PEXSI_INCLUDE_DIRS CP2K_PEXSI_FOUND)
