#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

file(READ "${PROJECT_SOURCE_DIR}/cmake/CompilerConfiguration.json"
     CP2K_COMPILER_CONFIGURATION NEWLINE_CONSUME)
foreach(
  __var
  CMAKE_Fortran_FLAGS
  CMAKE_Fortran_FLAGS_RELEASE
  CMAKE_Fortran_FLAGS_COVERAGE
  CMAKE_Fortran_FLAGS_DEBUG
  CMAKE_Fortran_FLAGS_DEBUG
  CMAKE_C_FLAGS
  CMAKE_C_FLAGS_RELEASE
  CMAKE_C_FLAGS_COVERAGE
  CMAKE_C_FLAGS_DEBUG)
  string(
    JSON
    cp2k_reading_value
    ERROR_VARIABLE
    cp2k_JSON_ERROR
    GET
    ${CP2K_COMPILER_CONFIGURATION}
    "compiler_id"
    "${CMAKE_Fortran_COMPILER_ID}"
    "${__var}")
  if(NOT cp2k_JSON_ERROR)
    set(${__var} "${${__var}} ${cp2k_reading_value}")
  endif()
endforeach()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    message(WARNING ${CMAKE_Fortran_FLAGS})
  endif()
endif()
#
