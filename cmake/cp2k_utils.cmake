
function(cp2k_set_default_paths _varname)
  #find_library should work when ${PACKAGE_ROOT} is given to cmake
  #(-DPACKAGE_ROOT=bla) but I use only one variable syntax CP2K_PACKAGE_PREFIX

  if (DEFINED ${_varname}_ROOT)
    set(CP2K_${_varname}_PREFIX
      "${${_varname}_ROOT}"
      PARENT_SCOPE)
  endif()

  # search common environment variables names
  if (NOT CP2K_${_varname}_PREFIX)
    foreach(__var ${_varname}_ROOT CRAY_${_varname}_ROOT OLCF_${_varname}_ROOT ${_varname}_PREFIX ${_varname}ROOT EB${_varname}ROOT)
	    if(DEFINED ENV{${__var}})
	      set(CP2K_${_varname}_PREFIX
          $ENV{${__var}}
          PARENT_SCOPE)
      endif()
    endforeach()

    # search for the default path
    if (NOT CP2K_${_varname}_PREFIX)
	    set(CP2K_${_varname}_PREFIX "/usr"  PARENT_SCOPE)
    endif()
  endif()
  mark_as_advanced(CP2K_${_varname}_PREFIX)
endfunction()

function(cp2k_find_libraries _package_name _library_name)
  string(TOUPPER ${_library_name} _library_name_upper)
  find_library(
    CP2K_${_package_name}_LIBRARIES
    NAMES "${_library_name}"
    PATHS "${CP2K_${_package_name}_PREFIX}"
    HINTS "${CP2K_${_package_name}_PREFIX}")
  if (CP2K_${_package_name}_LIBRARIES)
    set(CP2K_${_package_name}_LINK_LIBRARIES "${CP2K_${_package_name}_LIBRARIES}" PARENT_SCOPE)
    set(CP2K_${_package_name}_FOUND ON PARENT_SCOPE)
  endif()
  mark_as_advanced(CP2K_${_package_name}_LINK_LIBRARIES CP2K_${_package_name}_LIBRARIES)
endfunction()

function(cp2k_include_dirs _package_name _library_include_file)
  string(TOUPPER ${_package_name} _library_name_upper)
  find_path(
    CP2K_${_library_name_upper}_INCLUDE_DIRS
    NAMES "${_lib_include_file}"
    PATHS "${CP2K_${_package_name}_PREFIX}"
    HINTS "${CP2K_${_package_name}_PREFIX}"
    PATH_SUFFIXES "include" "include/${_pacakge_name}" "${_package_name}"
    NO_DEFAULT_PATH)
  mark_as_advanced(CP2K_${_library_name_upper}_INCLUDE_DIRS)
endfunction()
