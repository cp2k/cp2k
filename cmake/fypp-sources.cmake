#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

if(NOT TARGET fypp)
  add_custom_target(fypp) # common target for all fypp calls
endif()

# Use a system-provided fypp if available, otherwise the bundled one
find_program(
  FYPP_EXECUTABLE fypp
  DOC "The FYPP preprocessor"
  PATHS ../tools/build_utils)
if(NOT FYPP_EXECUTABLE)
  message(FATAL_ERROR "Failed to find the FYPP preprocessor.")
else()
  message(STATUS "FYPP preprocessor found.")
endif()

# https://gitlab.kitware.com/cmake/cmake/issues/18188
if((CMAKE_Fortran_COMPILER_ID STREQUAL "GNU") AND (CMAKE_VERSION
                                                   VERSION_GREATER_EQUAL 3.16))
  set(fypp_flags --line-numbering --line-marker-format=gfortran5)
endif()

function(ADD_FYPP_SOURCES OUTVAR)
  set(outfiles)

  foreach(f ${ARGN})
    # first we might need to make the input file absolute
    get_filename_component(f "${f}" ABSOLUTE)
    get_filename_component(ext "${f}" EXT)
    # get the relative path of the file to the current source dir
    file(RELATIVE_PATH rf "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
    # set the output filename of fypped sources
    set(of "${CMAKE_CURRENT_BINARY_DIR}/${rf}")

    # create the output directory if it doesn't exist
    get_filename_component(d "${of}" PATH)
    if(NOT IS_DIRECTORY "${d}")
      file(MAKE_DIRECTORY "${d}")
    endif()

    if("${f}" MATCHES ".F$")
      # append the output file to the list of outputs
      list(APPEND outfiles "${of}")
      # now add the custom command to generate the output file
      add_custom_command(
        OUTPUT "${of}"
        COMMAND ${Python_EXECUTABLE} ${FYPP_EXECUTABLE} ARGS ${fypp_flags}
                "${f}" "${of}"
        MAIN_DEPENDENCY "${f}"
        VERBATIM)
    elseif("${f}" MATCHES ".h$")
      # append the output file to the list of outputs
      list(APPEND outfiles "${of}")
      # now add the custom command to generate the output file
      add_custom_command(
        OUTPUT "${of}"
        COMMAND ${Python_EXECUTABLE} ${FYPP_EXECUTABLE} ARGS "-F" "${f}" "${of}"
        DEPENDS "${f}")
    else()
      configure_file("${f}" "${of}" COPYONLY)
    endif()
  endforeach()

  # build a custom target to fypp seperately (required for example by the doc
  # target) cmake 3.27.5 seems to have issues with these two next commands. I do
  # not use the fypp target anyway

  # add_custom_target("fypp_${OUTVAR}" DEPENDS ${outfiles})
  # add_dependencies(fypp "fypp_${OUTVAR}")

  # set the output list in the calling scope
  set(${OUTVAR}
      ${outfiles}
      PARENT_SCOPE)
endfunction()
