add_custom_target(fypp) # common target for all fypp calls

# Use a system-provided fypp if available, otherwise the bundled one
find_program(
  FYPP_EXECUTABLE fypp
  DOC "The FYPP preprocessor"
  PATHS ../tools/build_utils/fypp/bin)
if (NOT FYPP_EXECUTABLE)
  message(FATAL_ERROR "Failed to find the FYPP preprocessor.")
else ()
  message(STATUS "FYPP preprocessor found.")
endif ()

if ((CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    AND (CMAKE_GENERATOR STREQUAL "Ninja")
    AND (CMAKE_VERSION VERSION_GREATER_EQUAL 3.16))
  set(fypp_flags --line-numbering --line-marker-format=gfortran5)
elseif (CMAKE_BUILD_TYPE MATCHES COVERAGE)
  message(
    WARNING
      "Coverage build requested but your environment does not support Line Control directives in Fypp"
  )
  message(
    WARNING
      "You need CMake 3.16+, Ninja (CMake-patched) and gfortran 5+ for this to work!"
  )
  # otherwise the referenced lines in the Coverage report point to either the
  # original (unexpanded files) or to the Fypped sources which may then not be
  # picked up by the postprocessing tools. CMake 3.16+ and Ninja are needed
  # since older CMake (or CMake with make) are unable to parse Line Control
  # directives within line-continued USE stmts, see
  # https://gitlab.kitware.com/cmake/cmake/issues/18188
endif ()

function (ADD_FYPP_SOURCES OUTVAR)
  set(outfiles)

  foreach (f ${ARGN})
    # first we might need to make the input file absolute
    get_filename_component(f "${f}" ABSOLUTE)
    get_filename_component(ext "${f}" EXT)
    # get the relative path of the file to the current source dir
    file(RELATIVE_PATH rf "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
    # set the output filename of fypped sources
    set(of "${CMAKE_CURRENT_BINARY_DIR}/${rf}")

    # create the output directory if it doesn't exist
    get_filename_component(d "${of}" PATH)
    if (NOT IS_DIRECTORY "${d}")
      file(MAKE_DIRECTORY "${d}")
    endif ()

    if ("${f}" MATCHES ".F$")
      # append the output file to the list of outputs
      list(APPEND outfiles "${of}")
      # now add the custom command to generate the output file
      add_custom_command(
        OUTPUT "${of}"
        COMMAND ${Python_EXECUTABLE} ${FYPP_EXECUTABLE} ARGS ${fypp_flags}
                "${f}" "${of}"
        MAIN_DEPENDENCY "${f}"
        VERBATIM)
    elseif ("${f}" MATCHES ".h$")
      # append the output file to the list of outputs
      list(APPEND outfiles "${of}")
      # now add the custom command to generate the output file
      add_custom_command(
        OUTPUT "${of}"
        COMMAND ${Python_EXECUTABLE} ${FYPP_EXECUTABLE} ARGS "-F" "${f}" "${of}"
        DEPENDS "${f}")
    else ()
      configure_file("${f}" "${of}" COPYONLY)
    endif ()
  endforeach ()

  # build a custom target to fypp seperately (required for example by the doc
  # target)
  add_custom_target("fypp_${OUTVAR}" DEPENDS ${outfiles})
  add_dependencies(fypp "fypp_${OUTVAR}")

  # set the output list in the calling scope
  set(${OUTVAR}
      ${outfiles}
      PARENT_SCOPE)
endfunction ()
