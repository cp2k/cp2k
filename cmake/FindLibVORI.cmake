include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(LIBVORI)
cp2k_find_libraries(LIBVORI vori)
#cp2k_include_dirs(LIBVORI )

if(CP2K_LIBVORI_INCLUDE_DIRS)
  find_package_handle_standard_args(LibVORI DEFAULT_MSG CP2K_LIBVORI_LINK_LIBRARIES
    CP2K_LIBVORI_INCLUDE_DIRS)
else()
  find_package_handle_standard_args(LibVORI DEFAULT_MSG CP2K_LIBVORI_LINK_LIBRARIES)
endif()

if(NOT TARGET CP2K_VORI::vori)
  add_library(CP2K_VORI::vori INTERFACE IMPORTED)
  set_target_properties(CP2K_VORI::vori PROPERTIES INTERFACE_LINK_LIBRARIES
                                              "${CP2K_LIBVORI_LINK_LIBRARIES}")
  if(CP2K_LIBVORI_INCLUDE_DIRS)
    set_target_properties(CP2K_VORI::vori PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                                "${CP2K_LIBVORI_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBVORI_PREFIX CP2K_LIBVORI_INCLUDE_DIRS CP2K_LIBVORI_LINK_LIBRARIES CP2K_LIBVORI_LIBRARIES)
