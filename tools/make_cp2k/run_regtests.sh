#!/usr/bin/env bash
#
# run_regtests.sh: Run the CP2K regression tests (step 4 of make_cp2k.sh)
#
# Usage: run_regtests.sh CONFIG_FILE...
#
# The settings are read from the configuration files options.conf, deps.conf,
# and build.conf written by make_cp2k.sh and build_deps.sh.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"
activate_spack_env

# Launch regression test run
echo -e "\n*** Launching regression test run using the script ${INSTALL_PREFIX}/bin/run_tests\n"
if [[ "${TEST_COVERAGE}" == "yes" ]]; then
  # Print only a warning when the regression test is failing and continue with coverage analysis
  if ! ${LAUNCH_SCRIPT} run_tests; then
    echo -e "\nWARNING: The regression test run failed, but the coverage analysis will still be performed\n"
  fi
else
  ${LAUNCH_SCRIPT} run_tests
  EXIT_CODE=$?
  if ((EXIT_CODE != 0)); then
    echo -e "\nERROR: The regression test run failed with the error code ${EXIT_CODE}\n"
    exit "${EXIT_CODE}"
  fi
fi
