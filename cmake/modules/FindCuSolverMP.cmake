#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Copyright (c) 2022- ETH Zurich
#
# authors : Mathieu Taillefumier

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

find_package(ucc REQUIRED)

# First, find CuSolverMP library and headers
cp2k_set_default_paths(CUSOLVER_MP "CUSOLVER_MP")
cp2k_find_libraries(CUSOLVER_MP "cusolverMp")
cp2k_include_dirs(CUSOLVER_MP "cusolverMp.h")

# CuSolverMP 0.7+ uses NCCL for communication, older versions use Cal. We need
# to detect the version to require the correct communication library.
if(CP2K_CUSOLVER_MP_INCLUDE_DIRS)
  file(STRINGS "${CP2K_CUSOLVER_MP_INCLUDE_DIRS}/cusolverMp.h" _ver_major_line
       REGEX "^#define CUSOLVERMP_VER_MAJOR")
  file(STRINGS "${CP2K_CUSOLVER_MP_INCLUDE_DIRS}/cusolverMp.h" _ver_minor_line
       REGEX "^#define CUSOLVERMP_VER_MINOR")

  string(REGEX MATCH "[0-9]+" _ver_major "${_ver_major_line}")
  string(REGEX MATCH "[0-9]+" _ver_minor "${_ver_minor_line}")

  if(_ver_major STREQUAL "" OR _ver_minor STREQUAL "")
    message(FATAL_ERROR "Could not determine CuSolverMP version from header")
  endif()

  message(STATUS "Found CuSolverMP version: ${_ver_major}.${_ver_minor}")

  # CuSolverMP 0.7+ uses NCCL, older versions use Cal
  if(_ver_major GREATER 0 OR _ver_minor GREATER_EQUAL 7)
    find_package(Nccl REQUIRED)
    set(CP2K_CUSOLVERMP_USE_NCCL
        ON
        CACHE BOOL "CuSolverMP uses NCCL for communication")
  else()
    find_package(Cal REQUIRED)
    set(CP2K_CUSOLVERMP_USE_NCCL
        OFF
        CACHE BOOL "CuSolverMP uses Cal for communication")
  endif()
endif()

find_package_handle_standard_args(
  CuSolverMP DEFAULT_MSG CP2K_CUSOLVER_MP_LINK_LIBRARIES
  CP2K_CUSOLVER_MP_INCLUDE_DIRS)

if(CP2K_CUSOLVER_MP_FOUND AND NOT TARGET cp2k::CUSOLVER_MP::cusolver_mp)
  add_library(cp2k::CUSOLVER_MP::cusolver_mp INTERFACE IMPORTED)

  if(CP2K_CUSOLVERMP_USE_NCCL)
    set(_comm_lib "cp2k::NCCL::nccl")
  else()
    set(_comm_lib "cp2k::CAL::cal")
  endif()

  set_target_properties(
    cp2k::CUSOLVER_MP::cusolver_mp
    PROPERTIES INTERFACE_LINK_LIBRARIES
               "${CP2K_CUSOLVER_MP_LINK_LIBRARIES};${_comm_lib};cp2k::UCC::ucc")
  set_target_properties(
    cp2k::CUSOLVER_MP::cusolver_mp
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_CUSOLVER_MP_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "CuSolverMP requested, but not found")
endif()

mark_as_advanced(CP2K_CUSOLVER_MP_LINK_LIBRARIES)
mark_as_advanced(CP2K_CUSOLVER_MP_INCLUDE_DIRS)
mark_as_advanced(CP2K_CUSOLVER_MP_FOUND)
