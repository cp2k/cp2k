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

find_package(PkgConfig)

cp2k_set_default_paths(ARMPL "Armpl")

foreach(_var armpl_ilp64 armpl_lp64 armpl_ilp64_mp armpl_lp64_mp)
  string(TOUPPER ${_var} _var_up)
  cp2k_find_libraries(${_var_up} ${_var})
endforeach()

cp2k_include_dirs(ARMPL "armpl.h")

# Check for 64bit Integer support
if(CP2K_BLAS_INTERFACE MATCHES "64bits")
  set(CP2K_BLAS_armpl_LIB "armpl_ilp64")
else()
  set(CP2K_BLAS_armpl_LIB "armpl_lp64")
endif()

# Check for OpenMP support, VIA BLAS_VENDOR of Arm_mp or Arm_ipl64_mp
if(CP2K_BLAS_THREADING MATCHES "openmp")
  string(APPEND CP2K_BLAS_armpl_LIB "_mp")
endif()

# check if found
find_package_handle_standard_args(
  Armpl
  REQUIRED_VARS
    CP2K_ARMPL_INCLUDE_DIRS CP2K_ARMPL_LP64_LIBRARIES
    CP2K_ARMPL_LP64_MP_LIBRARIES CP2K_ARMPL_ILP64_LIBRARIES
    CP2K_ARMPL_ILP64_MP_LIBRARIES)

# add target to link against
if(CP2K_ARMPL_LP64_FOUND)

  if(NOT TARGET Armpl::armpl)
    add_library(CP2K::BLAS::Armpl::armpl INTERFACE IMPORTED)
    # now define an alias to the target library
    add_library(CP2K::BLAS::Armpl::blas ALIAS CP2K::BLAS::Armpl::armpl)
  endif()

  # we need to iniitialize the targets of each individual libraries only once.
  if(NOT TARGET CP2K::BLAS::Armpl::${_var})
    foreach(_var armpl_ilp64 armpl_lp64 armpl_ilp64_mp armpl_lp64_mp)
      string(TOUPPER "CP2K_${_var}_LINK_LIBRARIES" _var_up)
      add_library(CP2K::BLAS::Armpl::${_var} INTERFACE IMPORTED)
      set_property(
        TARGET CP2K::BLAS::Armpl::${_var} PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                   ${CP2K_ARMPL_INCLUDE_DIRS})
      set_property(TARGET CP2K::BLAS::Armpl::${_var}
                   PROPERTY INTERFACE_LINK_LIBRARIES "${${_var_up}}")
      endif()
    endforeach()
  endif()

  set_property(
    TARGET CP2K::BLAS::Armpl::armpl PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                             ${CP2K_ARMPL_INCLUDE_DIRS})
  set_property(
    TARGET CP2K::BLAS::Armpl::armpl
    PROPERTY INTERFACE_LINK_LIBRARIES
             "CP2K_${CP2K_BLAS_armpl_LIB}_LINK_LIBRARIES")
endif()

mark_as_advanced(CP2K_ARMPL_FOUND CP2K_ARMPL_LINK_LIBRARIES
                 CP2K_ARMPL_INCLUDE_DIRS)
