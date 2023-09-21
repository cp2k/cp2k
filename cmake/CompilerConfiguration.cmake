#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

file(READ "${PROJECT_SOURCE_DIR}/cmake/CompilerConfiguration.json"
     CP2K_COMPILER_CONFIGURATION NEWLINE_CONSUME)

# set(cp2k_supported_compiler_list "GNU|Intel|PGI|CRAY|NAG|Clang|AppleClang")
set(cp2k_supported_compiler_list "")
string(
  JSON
  cp2k_supported_compiler_list
  ERROR_VARIABLE
  cp2k_JSON_ERROR
  GET
  ${CP2K_COMPILER_CONFIGURATION}
  "compiler_id"
  "compiler_list")

function(cp2k_initialize_compiler_options _compiler_language _json_string)
  if("${CMAKE_${_compiler_language}_COMPILER_ID}" MATCHES
     ${cp2k_supported_compiler_list})
    foreach(
      __var
      CMAKE_${_compiler_language}_FLAGS
      CMAKE_${_compiler_language}_FLAGS_RELEASE
      CMAKE_${_compiler_language}_FLAGS_COVERAGE
      CMAKE_${_compiler_language}_FLAGS_DEBUG)

      string(
        JSON
        cp2k_reading_value
        ERROR_VARIABLE
        cp2k_JSON_ERROR
        GET
        ${_json_string}
        "compiler_id"
        "${CMAKE_${_compiler_language}_COMPILER_ID}"
        "${__var}")

      if(NOT cp2k_JSON_ERROR)
        set(${__var}
            "${${__var}} ${cp2k_reading_value}"
            PARENT_SCOPE)
      endif()
    endforeach()
  else()

    string(
      JSON
      cp2k_compiler_warning
      GET
      ${_json_string}
      "compiler_id"
      "warning"
      "${_compiler_language}")

    message(WARNING ${cp2k_compiler_warning})
    message("-- CMAKE_${_compiler_language}_COMPILER_ID: "
            ${CMAKE_${_compiler_language}_COMPILER_ID})
    message("-- CMAKE_${_compiler_language}_COMPILER full path: "
            ${CMAKE_${_compiler_language}_COMPILER})
  endif()
endfunction()

cp2k_initialize_compiler_options("Fortran" ${CP2K_COMPILER_CONFIGURATION})
cp2k_initialize_compiler_options("C" ${CP2K_COMPILER_CONFIGURATION})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # gcc 11 and above does not include -fallow-argument-mismatch
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()

  # the sanitizer reports leaks when cp2k is compiled with openmpi
  if((NOT (USE_MPI)) OR (NOT ("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES
                              "Open MPI")))
    set(CMAKE_Fortran_FLAGS_COVERAGE
        "${CMAKE_Fortran_FLAGS_COVERAGE} -fsanitize=leak")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG} -fsanitize=leak")
  endif()
endif()
