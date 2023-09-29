#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(FFTW3 "Fftw")
# Check if we can use PkgConfig
find_package(PkgConfig)

# First try with pkg
if(PKG_CONFIG_FOUND)
  pkg_search_module(CP2K_FFTW3 IMPORTED_TARGET GLOBAL fftw3)
  pkg_search_module(CP2K_FFTW3F IMPORTED_TARGET GLOBAL fftw3f)
  pkg_search_module(CP2K_FFTW3L IMPORTED_TARGET GLOBAL fftw3l)
  pkg_search_module(CP2K_FFTW3Q IMPORTED_TARGET GLOBAL fftw3q)
endif()

foreach(_lib fftw3 fftw3f fftw3l fftw3q)
  string(TOUPPER "${_lib}" __lib_up)
  if(NOT CP2K_${__lib_up}_FOUND)
    if(NOT ${_lib} MATCHES "fftw3")
      set(CP2K_${__lib_up}_ROOT "${CP2K_FFTW3_ROOT}")
    endif()
    cp2k_find_libraries("${__lib_up}" "${_lib}")
    if(NOT ${_lib} MATCHES "fftw3")
      unset(CP2K_${__lib_up}_ROOT CACHE)
    endif()
  endif()

  # OMP variant
  foreach(_subtype "mpi" "omp" "threads")
    string(TOUPPER "${_lib}_${_subtype}" _sub_lib)
    set(CP2K_${_sub_lib}_ROOT "${CP2K_FFTW3_ROOT}")
    cp2k_find_libraries("${_sub_lib}" "${_lib}_${_subtype}")
    unset(CP2K_${_sub_lib}_ROOT CACHE)
  endforeach()
endforeach()

if(NOT CP2K_FFTW3_INCLUDE_DIRS)
  cp2k_include_dirs(FFTW3 "fftw3.h;fftw3/fftw3.h")
endif()

if(CP2K_FFTW3_INCLUDE_DIRS)
  find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_FFTW3_INCLUDE_DIRS
                                    CP2K_FFTW3_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_FFTW3_LINK_LIBRARIES)
endif()

foreach(lib_name "fftw3" "fftw3l" "fftw3q" "fftw3f")
  string(TOUPPER "${lib_name}" __lib_name_up)

  if(CP2K_${__lib_name_up}_FOUND)
    if(NOT TARGET cp2k::FFTW3::${lib_name})
      add_library(cp2k::FFTW3::${lib_name} INTERFACE IMPORTED)
    endif()
    # we do not recheck if the libraries are found when pkg_config is
    # successful.
    set_target_properties(
      cp2k::FFTW3::${lib_name}
      PROPERTIES INTERFACE_LINK_LIBRARIES
                 "${CP2K_${__lib_name_up}_LINK_LIBRARIES}")

    if(CP2K_FFTW3_INCLUDE_DIRS)
      set_target_properties(
        cp2k::FFTW3::${lib_name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                            "${CP2K_FFTW3_INCLUDE_DIRS}")
    endif()

    foreach(sub_type "threads" "mpi" "omp")
      string(TOUPPER "${lib_name}_${sub_type}" __libs)
      if(CP2K_${__libs}_FOUND)
        if(NOT TARGET cp2k::FFTW3::${lib_name}_${sub_type})
          add_library(cp2k::FFTW3::${lib_name}_${sub_type} INTERFACE IMPORTED)
        endif()
        set_target_properties(
          cp2k::FFTW3::${lib_name}_${sub_type}
          PROPERTIES INTERFACE_LINK_LIBRARIES
                     "${CP2K_${__libs}_LINK_LIBRARIES}")
      endif()
    endforeach()
  endif()
endforeach()

set(CP2K_FFTW3_FOUND ON)
mark_as_advanced(
  CP2K_FFTW3_FOUND
  CP2K_FFTW3_ROOT
  CP2K_FFTW3_INCLUDE_DIRS
  CP2K_FFTW3_MPI
  CP2K_FFTW3_OMP
  CP2K_FFTW3_THREADS
  CP2K_FFTW3Q_OMP
  CP2K_FFTW3Q_THREADS
  CP2K_FFTW3F_MPI
  CP2K_FFTW3_OMP
  CP2K_FFTW3F_THREADS
  CP2K_FFTW3L_MPI
  CP2K_FFTW3L_OMP
  CP2K_FFTW3L_THREADS)
