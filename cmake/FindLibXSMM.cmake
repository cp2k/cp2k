# Copyright ETH Zurich 2022 -
#
include(FindPackageHandleStandardArgs)
include(cp2k_utils)
find_package(PkgConfig REQUIRED)

cp2k_set_default_paths(LIBXSMM)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(LIBXSMM IMPORTED_TARGET GLOBAL libxsmm)
  pkg_check_modules(LIBXSMMEXT IMPORTED_TARGET GLOBAL libxsmmext)
  pkg_check_modules(LIBXSMMF IMPORTED_TARGET GLOBAL libxsmmf)
  pkg_check_modules(LIBXSMMNOBLAS IMPORTED_TARGET GLOBAL libxsmmnoblas)
endif()

foreach (__lib libxsmm libxsmmf libxsmmext libxsmmnoblas)
  string(TOUPPER "${__lib}" __lib_search_up)
  if(NOT ${__lib_search_up}_FOUND)
    cp2k_find_libraries(${__lib_search_up} ${__lib})
  else()
    set(CP2K_${__lib_search_up}_LINK_LIBRARIES "${${__lib_search_up}_LINK_LIBRARIES}")
    set(CP2K_${__lib_search_up}_FOUND ON)
    set(CP2K_LIBXSMM_INCLUDE_DIRS "${LIBXSMM_INCLUDE_DIRS}")
  endif()
endforeach()

cp2k_include_dirs(LIBXSMM "libxsmm.h;libxsmm.mod;libxsmm/libxsmm.h")

if(CP2K_LIBXSMM_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibXSMM
    DEFAULT_MSG
    CP2K_LIBXSMM_INCLUDE_DIRS
    CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
    CP2K_LIBXSMMEXT_LINK_LIBRARIES
    CP2K_LIBXSMMF_LINK_LIBRARIES
    CP2K_LIBXSMM_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(
    LibXSMM DEFAULT_MSG CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES CP2K_LIBXSMMEXT_LINK_LIBRARIES 
    CP2K_LIBXSMMF_LINK_LIBRARIES CP2K_LIBXSMM_LINK_LIBRARIES)
  set(CP2K_LIBXSMM_INCLUDE_DIRS "/usr/include")
endif()

if(NOT TARGET CP2K_LibXSMM::libxsmm)
  foreach (__lib libxsmm libxsmmf libxsmmext libxsmmnoblas)
    string(TOUPPER "CP2K_${__lib}" __lib_search_up)
    if (${__lib_search_up}_FOUND AND NOT TARGET CP2K_LibXSMM::${__lib})
      add_library(CP2K_LibXSMM::${__lib} INTERFACE IMPORTED)
      set_target_properties(CP2K_LibXSMM::${__lib} PROPERTIES
        INTERFACE_LINK_LIBRARIES
        "${${__lib_search_up}_LINK_LIBRARIES}")
    endif()
  endforeach()
  if(CP2K_LIBXSMM_INCLUDE_DIRS)
    set_target_properties(
      CP2K_LibXSMM::libxsmm PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
      "${CP2K_LIBXSMM_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(CP2K_LIBXSMM_INCLUDE_DIRS
  CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
  CP2K_LIBXSMMEXT_LINK_LIBRARIES
  CP2K_LIBXSMMF_LINK_LIBRARIES
  CP2K_LIBXSMM_LINK_LIBRARIES)
