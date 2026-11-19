#!/usr/bin/env bash
#
# check_conventions.sh: Check the compliance with the CP2K coding conventions (step 4 of make_cp2k.sh)
#
# Usage: check_conventions.sh CONFIG_FILE...
#
# The settings are read from the configuration files options.conf, deps.conf,
# and build.conf written by make_cp2k.sh and build_deps.sh.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"
activate_spack_env

# Collect information from coding convention checks
"${CP2K_ROOT}"/tools/conventions/analyze_gfortran_ast.py "${CMAKE_BUILD_PATH}"/*.ast --jobs "${NUM_PROCS}" \
  &> "${CMAKE_BUILD_PATH}"/ast.issues
((VERBOSE > 0)) && cat "${CMAKE_BUILD_PATH}"/ast.issues
"${CP2K_ROOT}"/tools/conventions/analyze_gfortran_warnings.py "${CMAKE_BUILD_PATH}"/*.warn --jobs "${NUM_PROCS}" \
  &> "${CMAKE_BUILD_PATH}"/warn.issues
((VERBOSE > 0)) && cat "${CMAKE_BUILD_PATH}"/warn.issues
"${CP2K_ROOT}"/tools/conventions/summarize_issues.py --suppressions="${CP2K_ROOT}/tools/conventions/conventions.supp" "${CMAKE_BUILD_PATH}"/*.issues
cat << *** > "${INSTALL_PREFIX}"/bin/summarize_issues
#!/bin/bash
${CP2K_ROOT}/tools/conventions/summarize_issues.py --suppressions=${CP2K_ROOT}/tools/conventions/conventions.supp ${CMAKE_BUILD_PATH}/*.issues
***
chmod 750 "${INSTALL_PREFIX}"/bin/summarize_issues
