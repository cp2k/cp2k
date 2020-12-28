# =================================================================================================
# LCOV - COVERAGE REPORTS GENERATION
find_program(
  LCOV_EXE lcov
  DOC "path to the lcov executable (required to generate coverage reports)")

find_program(
  GENHTML_EXE genhtml
  DOC "path to the genhtml executable (required to generate HTML coverage reports)"
)

add_custom_target(
  cov-info
  COMMAND
    "${LCOV_EXE}" --directory "${CMAKE_BINARY_DIR}" --base-dir
    "${CMAKE_SOURCE_DIR}" --no-external --capture --output-file coverage.info
  COMMAND "${LCOV_EXE}" --list coverage.info
  VERBATIM
  BYPRODUCTS coverage.info
  COMMENT "Generate coverage.info using lcov")

add_custom_target(
  cov-genhtml
  COMMAND "${GENHTML_EXE}" coverage.info --output-directory cov-html
  COMMENT
    "Generate a HTML-based coverage report using lcov in ${CMAKE_BINARY_DIR}/cov-html"
  VERBATIM) # Note: this directory will not be cleaned by `cmake --build .. --
            # clean`
add_dependencies(cov-genhtml cov-info)
