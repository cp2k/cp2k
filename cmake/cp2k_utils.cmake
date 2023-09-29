#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

function(cp2k_set_default_paths _varname _package_name)
  # find_library should work when ${PACKAGE_ROOT} is given to cmake
  # (-DPACKAGE_ROOT=bla) but I use only one variable syntax CP2K_PACKAGE_PREFIX
  set(CP2K_${_varname}_PREFIX_TMP "")
  if(DEFINED ${_package_name}_ROOT)
    set(CP2K_${_varname}_PREFIX_TMP "${${_varname}_ROOT}")
  endif()

  # search common environment variables names
  if(NOT CP2K_${_varname}_PREFIX_TMP)
    foreach(
      __var
      ${_varname}_ROOT
      CRAY_${_varname}_PREFIX_DIR
      CRAY_${_varname}_ROOT
      OLCF_${_varname}_ROOT
      ${_varname}_PREFIX
      ${_varname}ROOT
      ${_varname}_HOME
      EB${_varname}ROOT)
      if(DEFINED ENV{${__var}})
        set(CP2K_${_varname}_PREFIX_TMP $ENV{${__var}})
      endif()
    endforeach()

    # search for the default path
    if(NOT CP2K_${_varname}_PREFIX_TMP)
      set(CP2K_${_varname}_PREFIX_TMP "/usr")
    endif()
  endif()
  set(CP2K_${_varname}_ROOT
      "${CP2K_${_varname}_PREFIX_TMP}"
      PARENT_SCOPE)

  unset(CP2K_${_varname}_PREFIX_TMP CACHE)
endfunction()

function(cp2k_find_libraries _package_name _library_name)
  string(TOUPPER ${_library_name} _library_name_upper)

  find_library(
    CP2K_${_package_name}_LIBRARIES_TMP
    NAMES ${_library_name}
    PATHS "${CP2K_${_package_name}_ROOT}"
    PATH_SUFFIXES "lib" "lib64")
  if(CP2K_${_package_name}_LIBRARIES_TMP)
    set(CP2K_${_package_name}_LINK_LIBRARIES
        "${CP2K_${_package_name}_LIBRARIES_TMP}"
        PARENT_SCOPE)
    set(CP2K_${_package_name}_LIBRARIES
        "${CP2K_${_package_name}_LIBRARIES_TMP}"
        PARENT_SCOPE)
    set(CP2K_${_package_name}_FOUND
        ON
        PARENT_SCOPE)
  endif()

  unset(CP2K_${_package_name}_LIBRARIES_TMP CACHE)
endfunction()

function(cp2k_include_dirs _package_name _library_include_file)
  find_path(
    CP2K_${_package_name}_INCLUDE_DIRS_TMP
    NAMES ${_library_include_file}
    PATHS "${CP2K_${_package_name}_ROOT}"
    HINTS "${CP2K_${_package_name}_ROOT}"
    PATH_SUFFIXES "include" "include/${_pacakge_name}" "${_package_name}")

  set(CP2K_${_package_name}_INCLUDE_DIRS
      "${CP2K_${_package_name}_INCLUDE_DIRS_TMP}"
      PARENT_SCOPE)
  unset(CP2K_${_package_name}_INCLUDE_DIRS_TMP CACHE)
endfunction()

function(cp2k_compare_src_with_list _list_files _extension)
  file(
    GLOB_RECURSE _test_list
    RELATIVE "${CMAKE_SOURCE_DIR}/src"
    "${_extension}")
  # message(STATUS "ref : ${_list_files}") message(STATUS "search :
  # ${_test_list}")

  set(_test_list_1 ${_test_list})
  # remove all item common to the two lists
  list(REMOVE_ITEM _test_list_1 ${_list_files})
  list(REMOVE_ITEM _list_files ${_test_list})

  # now _list_files and _test_list_1 should both be empty.
  list(LENGTH _test_list_1 found_list_size_)
  list(LENGTH _list_files list_size_)

  if((list_size_ GREATER 0) OR (found_list_size_ GREATER 0))
    message(
      FATAL_ERROR
        "The files registered in CMakeLists.txt and the files found with\n"
        "the extension ${_extension} do not match.\n\n"
        "Your src directory likely contains files that were either\n"
        "renamed/added/deleted or forgotten to be deleted.\n"
        "This list of files (from the CMakeLists.txt) was either deleted or renamed\n\n"
        "${_list_files}\n\n"
        "while a direct search returned these additional files (to be added or removed)\n\n"
        "${_test_list_1}\n\n"
        "Either add these files to CMakeLists.txt or remove them.")
  endif()
  set(_test_list "")
  set(_test_list_1 "")
  set(list_size_ 0)
  set(found_list_size_ 0)
endfunction()

macro(cp2k_FindPackage name)
  #[===[.md
    # cp2k_FindPackage

    A compatibility macro that links `find_package(CONFIG)` packages with `pkg-config`. This should only
    be called within the `Find<PackageName>.cmake` file.

    Note: Version range syntax is not supported for pkg-config searching. Only the lower bound will be respected.
          Taken from Octopus repository.
    ]===]

  set(ARGS_Options "")
  set(ARGS_OneValue "")
  set(ARGS_MultiValue "NAMES;PKG_MODULE_NAMES;PKG_MODULE_SPECS")
  cmake_parse_arguments(ARGS "${ARGS_Options}" "${ARGS_OneValue}"
                        "${ARGS_MultiValue}" ${ARGN})

  # First try to find native <PackageName>Config.cmake Build the arguments
  # COMPONENTS
  set(_comp_args)
  set(_opt_comp_args)
  if(DEFINED ${name}_FIND_COMPONENTS)
    list(APPEND _comp_args COMPONENTS)
    foreach(_comp IN LISTS ${name}_FIND_COMPONENTS)
      if(${name}_FIND_REQUIRED_${_comp})
        list(APPEND _comp_args ${_comp})
      else()
        if(NOT DEFINED _opt_comp_args)
          list(APPEND _opt_comp_args OPTIONAL_COMPONENTS)
        endif()
        list(APPEND _opt_comp_args ${_comp})
      endif()
    endforeach()
  endif()

  # Version Try range format first, otherwise use the default
  set(_version_args ${${name}_FIND_VERSION_RANGE})
  if(NOT DEFINED _version_args)
    set(_version_args ${${name}_FIND_VERSION})
  endif()
  if(${name}_FIND_VERSION_EXACT)
    list(APPEND _version_args EXACT)
  endif()

  # QUIET
  set(_quiet_arg)
  if(${name}_FIND_QUIETLY)
    list(APPEND _quiet_arg QUIET)
  endif()

  # REQUIRED
  set(_required_arg)
  if(${name}_FIND_REQUIRED)
    list(APPEND _required_arg REQUIRED)
  endif()

  # REGISTRY_VIEW
  set(_registry_view_arg)
  if(${name}_FIND_REGISTRY_VIEW)
    list(APPEND _registry_view REGISTRY_VIEW ${${name}_FIND_REGISTRY_VIEW})
  endif()

  # NAMES
  set(_names_args)
  if(DEFINED ARGS_NAMES)
    list(APPEND _names_args NAMES ${ARGS_NAMES})
  endif()

  # Try <PackageName>Config.cmake
  find_package(
    ${name}
    ${_version_args}
    ${_quiet_arg}
    CONFIG
    ${_comp_args}
    ${_opt_comp_args}
    ${_registry_view_arg}
    ${_names_args})

  if(NOT ${name}_FOUND)
    # Try pkg-config next Construct the moduleSpec to search for
    if(NOT DEFINED ARGS_PKG_MODULE_SPECS)
      if(NOT DEFINED ARGS_PKG_MODULE_NAMES)
        set(ARGS_PKG_MODULE_NAMES ${name})
      endif()
      if(DEFINED ${name}_FIND_VERSION_RANGE)
        # Can only parse the minimum requirement
        foreach(_pkg_name IN LISTS ARGS_PKG_MODULE_NAMES)
          list(APPEND ARGS_PKG_MODULE_SPECS
               "${_pkg_name}>=${${name}_FIND_VERSION_MIN}")
        endforeach()
      elseif({${name}_FIND_VERSION_EXACT)
        # Requesting exact version
        foreach(_pkg_name IN LISTS ARGS_PKG_MODULE_NAMES)
          list(APPEND ARGS_PKG_MODULE_SPECS
               "${_pkg_name}=${${name}_FIND_VERSION}")
        endforeach()
      elseif(DEFINED ${name}_FIND_VERSION)
        # Otherwise treat the request as minimum requirement
        foreach(_pkg_name IN LISTS ARGS_PKG_MODULE_NAMES)
          list(APPEND ARGS_PKG_MODULE_SPECS
               "${_pkg_name}>=${${name}_FIND_VERSION}")
        endforeach()
      else()
        # Fallthrough if no version is required
        foreach(_pkg_name IN LISTS ARGS_PKG_MODULE_NAMES)
          list(APPEND ARGS_PKG_MODULE_SPECS "${_pkg_name}")
        endforeach()
      endif()
    endif()
    # Call pkg-config
    if(CMAKE_VERSION VERSION_LESS 3.28)
      # https://gitlab.kitware.com/cmake/cmake/-/issues/25228
      set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
    endif()
    if(CMAKE_VERSION VERSION_LESS 3.22)
      # Back-porting
      # https://gitlab.kitware.com/cmake/cmake/-/merge_requests/6345
      set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
    endif()
    pkg_search_module(${name} ${_required_arg} ${_quiet_arg} IMPORTED_TARGET
                      ${ARGS_PKG_MODULE_SPECS})
    # Mark the package as found by pkg-config
    if(${name}_FOUND)
      set(${name}_PKGCONFIG True)
    endif()
  endif()

  # Sanitize local variables in order to not contaminate future calls
  set(ARGS_Options)
  set(ARGS_OneValue)
  set(ARGS_MultiValue)
  set(ARGS_UNPARSED_ARGUMENTS)
  set(ARGS_NAMES)
  set(ARGS_PKG_MODULE_NAMES)
  set(ARGS_PKG_MODULE_SPECS)
  set(_version_args)
  set(_quiet_arg)
  set(_comp_args)
  set(_opt_comp_args)
  set(_registry_view_arg)
  set(_names_args)
  set(_pkg_name)
endmacro()
