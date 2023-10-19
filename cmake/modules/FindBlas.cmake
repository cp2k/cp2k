#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

if(NOT
   (CMAKE_C_COMPILER_LOADED
    OR CMAKE_CXX_COMPILER_LOADED
    OR CMAKE_Fortran_COMPILER_LOADED))
  message(FATAL_ERROR "FindBLAS requires Fortran, C, or C++ to be enabled.")
endif()

if(NOT CP2K_CONFIG_PACKAGE)
  set(CP2K_BLAS_VENDOR_LIST
      # cmake-format: sortable
      "auto"
      "MKL"
      "OpenBLAS"
      "SCI"
      "GenericBLAS"
      "Armpl"
      "FlexiBLAS"
      "Atlas"
      "NVHPCBlas"
      "CUSTOM")

  set(__BLAS_VENDOR_LIST ${CP2K_BLAS_VENDOR_LIST})
  list(REMOVE_ITEM __BLAS_VENDOR_LIST "auto")
  list(REMOVE_ITEM __BLAS_VENDOR_LIST "CUSTOM")

  # set(CP2K_BLAS_VENDOR "auto" CACHE STRING "Blas library for computations on
  # host")
  set_property(CACHE CP2K_BLAS_VENDOR PROPERTY STRINGS ${CP2K_BLAS_VENDOR_LIST})

  if(NOT ${CP2K_BLAS_VENDOR} IN_LIST CP2K_BLAS_VENDOR_LIST)
    message(FATAL_ERROR "Invalid Host BLAS backend")
  endif()

  set(CP2K_BLAS_THREAD_LIST
      # cmake-format: sortable
      "sequential" "thread" "gnu-thread" "intel-thread" "tbb-thread" "openmp")

  set(CP2K_BLAS_THREADING
      "sequential"
      CACHE STRING "threaded blas library")
  set_property(CACHE CP2K_BLAS_THREADING PROPERTY STRINGS
                                                  ${CP2K_BLAS_THREAD_LIST})

  if(NOT ${CP2K_BLAS_THREADING} IN_LIST CP2K_BLAS_THREAD_LIST)
    message(FATAL_ERROR "Invalid threaded BLAS backend")
  endif()

  set(CP2K_BLAS_INTERFACE_BITS_LIST "32bits" "64bits")
  set(CP2K_BLAS_INTERFACE
      "32bits"
      CACHE STRING
            "32 bits integers are used for indices, matrices and vectors sizes")
  set_property(CACHE CP2K_BLAS_INTERFACE
               PROPERTY STRINGS ${CP2K_BLAS_INTERFACE_BITS_LIST})

  if(NOT ${CP2K_BLAS_INTERFACE} IN_LIST CP2K_BLAS_INTERFACE_BITS_LIST)
    message(
      FATAL_ERROR
        "Invalid parameters. Blas and lapack can exist in two flavors 32 or 64 bits interfaces (relevant mostly for mkl)"
    )
  endif()

  set(CP2K_BLAS_FOUND FALSE)

  # first check for a specific implementation if requested

  if(NOT CP2K_BLAS_VENDOR MATCHES "auto|CUSTOM")
    if(DEFINED CP2K_BLAS_LINK_LIBRARIES)
      set(CP2K_BLAS_FOUND TRUE)
    else()
      find_package(${CP2K_BLAS_VENDOR} REQUIRED)
      if(TARGET cp2k::BLAS::${CP2K_BLAS_VENDOR}::blas)
        get_target_property(
          CP2K_BLAS_INCLUDE_DIRS cp2k::BLAS::${CP2K_BLAS_VENDOR}::blas
          INTERFACE_INCLUDE_DIRECTORIES)
        get_target_property(
          CP2K_BLAS_LINK_LIBRARIES cp2k::BLAS::${CP2K_BLAS_VENDOR}::blas
          INTERFACE_LINK_LIBRARIES)
        set(CP2K_BLAS_FOUND TRUE)
      endif()
    endif()
  else()
    if(CP2K_BLAS_VENDOR MATCHES "CUSTOM" AND NOT DEFINED
                                             CP2K_BLAS_LINK_LIBRARIES)
      message(
        FATAL_ERROR
          "Setting CP2K_BLAS_VENDOR=CUSTOM imply setting CP2K_BLAS_LINK_LIBRARIES\n and CP2K_LAPACK_LINK_LIBRARIES to the right libraries. See the README_cmake.md for more details"
      )
    endif()

    if(DEFINED CP2K_BLAS_LINK_LIBRARIES)
      set(CP2K_BLAS_FOUND TRUE)
    else()
      # search for any blas implementation and exit immediately if one is found.
      # we could also give a full list of found implementation and let the user
      # choose which implementation to use
      foreach(_libs ${__BLAS_VENDOR_LIST})
        # I exclude the first item of the list
        find_package(${_libs})
        if(TARGET cp2k::BLAS::${_libs}::blas)
          get_target_property(CP2K_BLAS_INCLUDE_DIRS cp2k::BLAS::${_libs}::blas
                              INTERFACE_INCLUDE_DIRECTORIES)
          get_target_property(
            CP2K_BLAS_LINK_LIBRARIES cp2k::BLAS::${_libs}::blas
            INTERFACE_LINK_LIBRARIES)
          set(CP2K_BLAS_VENDOR "${_libs}")
          set(CP2K_BLAS_FOUND TRUE)
          break()
        endif()
      endforeach()
    endif()
  endif()
else()
  set(CP2K_BLAS_FOUND ON)
endif()

# we exclude the CP2K_BLAS_INCLUDE_DIRS from the list of mandatory variables as
# having the fortran interface is usually enough. C, C++ and others languages
# might require this information though

find_package_handle_standard_args(
  Blas REQUIRED_VARS CP2K_BLAS_LINK_LIBRARIES CP2K_BLAS_VENDOR CP2K_BLAS_FOUND)

if(NOT TARGET cp2k::BLAS::blas)
  add_library(cp2k::BLAS::blas INTERFACE IMPORTED)
endif()

set_target_properties(cp2k::BLAS::blas PROPERTIES INTERFACE_LINK_LIBRARIES
                                                  "${CP2K_BLAS_LINK_LIBRARIES}")

if(CP2K_BLAS_INCLUDE_DIRS)
  set_target_properties(
    cp2k::BLAS::blas PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                "${CP2K_BLAS_INCLUDE_DIRS}")
endif()

mark_as_advanced(CP2K_BLAS_INCLUDE_DIRS)
mark_as_advanced(CP2K_BLAS_LINK_LIBRARIES)
mark_as_advanced(CP2K_BLAS_VENDOR)
mark_as_advanced(CP2K_BLAS_FOUND)
