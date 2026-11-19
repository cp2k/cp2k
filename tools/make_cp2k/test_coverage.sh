#!/usr/bin/env bash
#
# test_coverage.sh: Analyse the code coverage and generate a coverage report (step 4 of make_cp2k.sh)
#
# Usage: test_coverage.sh CONFIG_FILE...
#
# The settings are read from the configuration files options.conf, deps.conf,
# and build.conf written by make_cp2k.sh and build_deps.sh.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"
activate_spack_env

# Analyse code coverage of the regression test run and generate a coverage report
if [[ "${IN_CONTAINER}" == "yes" ]]; then
  COVERAGE_OUTPUT_DIR="/workspace/artifacts/coverage"
else
  COVERAGE_OUTPUT_DIR="${INSTALL_PREFIX}"/coverage
fi
mkdir -p "${COVERAGE_OUTPUT_DIR}"
COVERAGE_OUTPUT_FILE="${COVERAGE_OUTPUT_DIR}/coverage.info"
lcov --directory "${CMAKE_BUILD_PATH}/src" \
  --exclude "${SPACK_BUILD_PATH}/*" \
  --exclude "/usr/*" \
  --capture \
  --keep-going \
  --output-file "${COVERAGE_OUTPUT_FILE}" &> "${COVERAGE_OUTPUT_DIR}"/lcov.log
# Print coverage summary
lcov --summary "${COVERAGE_OUTPUT_FILE}"
genhtml --output-directory "${COVERAGE_OUTPUT_DIR}" --keep-going --title "CP2K Regtests (${CP2K_REVISION})" \
  "${COVERAGE_OUTPUT_FILE}" &> "${COVERAGE_OUTPUT_DIR}"/genhtml.log
# Create plot data
LINE_COV=$(lcov --summary "${COVERAGE_OUTPUT_FILE}" | grep lines | awk '{print substr($2, 1, length($2)-1)}')
FUNC_COV=$(lcov --summary "${COVERAGE_OUTPUT_FILE}" | grep funct | awk '{print substr($2, 1, length($2)-1)}')
echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
echo "PlotPoint: name='lines', plot='cov', label='Lines', y=${LINE_COV}, yerr=0"
echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=${FUNC_COV}, yerr=0"
