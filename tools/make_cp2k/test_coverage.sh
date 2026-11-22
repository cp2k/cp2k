#!/usr/bin/env bash
#
# test_coverage.sh: Analyse the code coverage and generate a coverage report
#
# Usage: test_coverage.sh [OPTIONS] [CONFIG_FILE...]
#
# This script is used for the coverage tests of both the Spack-based and the
# toolchain-based CP2K builds:
#
#  - With configuration files (options.conf, deps.conf, and build.conf written by
#    make_cp2k.sh and build_deps.sh), the defaults are taken from these files and
#    the Spack environment with the CP2K dependencies is activated (step 4 of
#    make_cp2k.sh, which runs the regression tests separately)
#  - Without configuration files, at least the flag --build_dir is required and
#    the environment can be set up with the flag --setup, e.g. for a toolchain build:
#
#      test_coverage.sh --setup /opt/cp2k-toolchain/install/setup \
#                       --regtest "--ompthreads=1 --keepalive /opt/cp2k/build/bin psmp" \
#                       --build_dir /opt/cp2k/build/src --exclude "/opt/cp2k-toolchain/*"
#
# The flags override the defaults from the configuration files.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

CP2K_ROOT_DIR="$(realpath "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/../..")"

# shellcheck source=/dev/null
source "${CP2K_ROOT_DIR}/tools/make_cp2k/common.sh"

print_usage() {
  echo "Usage: $(basename "${0}") [OPTIONS] [CONFIG_FILE...]"
  echo ""
  echo "Options:"
  echo " -b, --build_dir DIR    : Build folder with the coverage data (.gcda files)"
  echo "                          (default: \${CMAKE_BUILD_PATH}/src)"
  echo " -o, --output_dir DIR   : Folder for the coverage report"
  echo "                          (default: /workspace/artifacts/coverage if \${IN_CONTAINER} is yes,"
  echo "                          otherwise \${INSTALL_PREFIX}/coverage or ./coverage)"
  echo " -x, --exclude PATTERN  : Exclude the files matching PATTERN from the coverage analysis"
  echo "                          (can be repeated, /usr/* and \${SPACK_BUILD_PATH}/* are always excluded)"
  echo " -t, --title TITLE      : Title of the HTML coverage report"
  echo "                          (default: CP2K Regtests (<CP2K revision>))"
  echo " -s, --setup FILE       : Source FILE to set up the environment, e.g. the toolchain setup file"
  echo " -r, --regtest ARGS     : Run the regression tests with tests/do_regtest.py ARGS before"
  echo "                          the coverage analysis"
  echo " -gt, --gcov_timeout T  : Remove the coverage data files for which gcov takes longer than T"
  echo "                          (e.g. 10s) to bypass gcov hanging on some files (default: no check)"
  echo " -h, --help             : Print this help information"
}

# Parse the flags
BUILD_DIR=""
OUTPUT_DIR=""
TITLE=""
SETUP_FILE=""
REGTEST_ARGS=""
GCOV_TIMEOUT=""
EXCLUDES=("/usr/*")
CONFIG_FILES=()
while (($# > 0)); do
  case "${1}" in
    -b | --build_dir | -o | --output_dir | -x | --exclude | -t | --title | -s | --setup | \
      -r | --regtest | -gt | --gcov_timeout)
      if (($# < 2)) || [[ -z "${2}" ]]; then
        echo "ERROR: No argument found for flag \"${1}\""
        exit 1
      fi
      case "${1}" in
        -b | --build_dir) BUILD_DIR="${2}" ;;
        -o | --output_dir) OUTPUT_DIR="${2}" ;;
        -x | --exclude) EXCLUDES+=("${2}") ;;
        -t | --title) TITLE="${2}" ;;
        -s | --setup) SETUP_FILE="${2}" ;;
        -r | --regtest) REGTEST_ARGS="${2}" ;;
        -gt | --gcov_timeout) GCOV_TIMEOUT="${2}" ;;
      esac
      shift 2
      ;;
    -h | --help)
      print_usage
      exit 0
      ;;
    -*)
      echo "ERROR: Unknown flag \"${1}\" found"
      print_usage
      exit 1
      ;;
    *)
      CONFIG_FILES+=("${1}")
      shift 1
      ;;
  esac
done

# Set up the environment and take the defaults from the configuration files if provided
if [[ -n "${SETUP_FILE}" ]]; then
  if [[ ! -f "${SETUP_FILE}" ]]; then
    echo "ERROR: The setup file ${SETUP_FILE} was not found"
    exit 1
  fi
  # The setup file is not necessarily safe for unset variables
  set +u
  # shellcheck source=/dev/null
  source "${SETUP_FILE}"
  set -u
fi
if ((${#CONFIG_FILES[@]} > 0)); then
  load_config "${CONFIG_FILES[@]}"
  activate_spack_env
  [[ -n "${BUILD_DIR}" ]] || BUILD_DIR="${CMAKE_BUILD_PATH}/src"
  [[ -n "${SPACK_BUILD_PATH:-}" ]] && EXCLUDES+=("${SPACK_BUILD_PATH}/*")
elif [[ -z "${BUILD_DIR}" ]]; then
  echo "ERROR: No configuration file and no build folder (flag --build_dir) specified"
  print_usage
  exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
  if [[ "${IN_CONTAINER:-}" == "yes" ]]; then
    OUTPUT_DIR="/workspace/artifacts/coverage"
  elif [[ -n "${INSTALL_PREFIX:-}" ]]; then
    OUTPUT_DIR="${INSTALL_PREFIX}/coverage"
  else
    OUTPUT_DIR="${PWD}/coverage"
  fi
fi
if [[ -z "${TITLE}" ]]; then
  if [[ -z "${CP2K_REVISION:-}" ]]; then
    CP2K_REVISION="$("${CP2K_ROOT_DIR}"/tools/build_utils/get_revision_number "${CP2K_ROOT_DIR}"/src)"
  fi
  TITLE="CP2K Regtests (${CP2K_REVISION})"
fi

if [[ ! -d "${BUILD_DIR}" ]]; then
  echo "ERROR: The build folder ${BUILD_DIR} (flag --build_dir or CMAKE_BUILD_PATH/src) was not found"
  exit 1
fi
BUILD_DIR="$(realpath "${BUILD_DIR}")"
for tool in gcov lcov genhtml; do
  if ! command -v "${tool}" &> /dev/null; then
    echo "ERROR: The command ${tool} required for the coverage analysis was not found"
    exit 1
  fi
done

# Run the regression tests if requested
if [[ -n "${REGTEST_ARGS}" ]]; then
  # Disable core dumps as they can take a very long time to write
  ulimit -c 0
  # Check available shared memory needed for MPI inter-process communication
  SHM_AVAIL=$(df --output=avail -m /dev/shm 2> /dev/null | tail -1 | tr -dc '0-9')
  if ((${SHM_AVAIL:-0} < 1024)); then
    echo "ERROR: Only ${SHM_AVAIL:-0} MB of shared memory available in /dev/shm, but at least 1024 MB are needed"
    echo "       (use --shm-size=1g with docker or podman)"
    exit 1
  fi
  echo -e "\n========== Running Regtests =========="
  read -ra REGTEST_ARGV <<< "${REGTEST_ARGS}"
  # A failure of the regression tests does not prevent the coverage analysis
  (cd "${CP2K_ROOT_DIR}" && ./tests/do_regtest.py "${REGTEST_ARGV[@]}")
  EXIT_CODE=$?
  if ((EXIT_CODE != 0)); then
    echo "WARNING: The regression tests failed (exit code ${EXIT_CODE}), the coverage analysis is continued anyway"
  fi
fi

# Check for coverage data
if [[ -z "$(find "${BUILD_DIR}" -name "*.gcda" -print -quit)" ]]; then
  echo "ERROR: No coverage data (.gcda files) found in ${BUILD_DIR}."
  echo "       Was CP2K built with coverage instrumentation and were the tests run?"
  exit 1
fi

echo -e "\n========== Analysing Coverage =========="

# Remove the coverage data files for which gcov takes too long, since gcov gets
# stuck on some files (maybe related to https://bugs.launchpad.net/gcc-arm-embedded/+bug/1694644)
if [[ -n "${GCOV_TIMEOUT}" ]]; then
  GCOV_WORK_DIR="$(mktemp -d)"
  find "${BUILD_DIR}" -name "*.gcda" -print0 | while read -r -d $'\0' gcda_file; do
    (cd "${GCOV_WORK_DIR}" && timeout "${GCOV_TIMEOUT}" gcov "${gcda_file}" &> /dev/null)
    # The exit code of timeout is 124 if the time limit was exceeded
    if (($? == 124)); then
      echo "Skipping ${gcda_file} because gcov took longer than ${GCOV_TIMEOUT}"
      rm -f "${gcda_file}"
    fi
  done
  rm -rf "${GCOV_WORK_DIR}"
fi

# Analyse the code coverage of the regression test run and generate a coverage report
mkdir -p "${OUTPUT_DIR}"
COVERAGE_OUTPUT_FILE="${OUTPUT_DIR}/coverage.info"
rm -f "${COVERAGE_OUTPUT_FILE}"
LCOV_EXCLUDE_FLAGS=()
for pattern in "${EXCLUDES[@]}"; do
  LCOV_EXCLUDE_FLAGS+=(--exclude "${pattern}")
done
lcov --directory "${BUILD_DIR}" \
  "${LCOV_EXCLUDE_FLAGS[@]}" \
  --capture \
  --keep-going \
  --output-file "${COVERAGE_OUTPUT_FILE}" &> "${OUTPUT_DIR}"/lcov.log
# lcov may also return a non-zero exit code for minor issues like unused exclude patterns
EXIT_CODE=$?
if ! grep -q "^SF:" "${COVERAGE_OUTPUT_FILE}" 2> /dev/null; then
  echo "ERROR: The coverage data capture with lcov failed (see ${OUTPUT_DIR}/lcov.log)"
  exit 1
elif ((EXIT_CODE != 0)); then
  echo "WARNING: lcov reported errors (exit code ${EXIT_CODE}, see ${OUTPUT_DIR}/lcov.log)"
fi
# Print coverage summary
lcov --summary "${COVERAGE_OUTPUT_FILE}"
rm -f "${OUTPUT_DIR}"/index.html
genhtml --output-directory "${OUTPUT_DIR}" --keep-going --title "${TITLE}" \
  "${COVERAGE_OUTPUT_FILE}" &> "${OUTPUT_DIR}"/genhtml.log
EXIT_CODE=$?
if [[ ! -s "${OUTPUT_DIR}"/index.html ]]; then
  echo "ERROR: The generation of the HTML coverage report with genhtml failed (see ${OUTPUT_DIR}/genhtml.log)"
  exit 1
elif ((EXIT_CODE != 0)); then
  echo "WARNING: genhtml reported errors (exit code ${EXIT_CODE}, see ${OUTPUT_DIR}/genhtml.log)"
fi
# Create plot data
LINE_COV=$(lcov --summary "${COVERAGE_OUTPUT_FILE}" | grep lines | awk '{print substr($2, 1, length($2)-1)}')
FUNC_COV=$(lcov --summary "${COVERAGE_OUTPUT_FILE}" | grep funct | awk '{print substr($2, 1, length($2)-1)}')
echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
echo "PlotPoint: name='lines', plot='cov', label='Lines', y=${LINE_COV}, yerr=0"
echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=${FUNC_COV}, yerr=0"
