# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Optional regression tests, coverage report, and performance run
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

# Set image tag if available
export IMAGE_TAG=${IMAGE_TAG:-<IMAGE ID>}

# Optionally, launch test run
if [[ "${RUN_TEST}" == "yes" ]]; then
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
else
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    if ((CUDA_SM_CODE > 0)); then
      DEVICE_FLAG=" --device nvidia.com/gpu=all"
    else
      DEVICE_FLAG=""
    fi
    echo ""
    echo "*** A regression test run can be launched with"
    echo "    podman run -it${DEVICE_FLAG} --rm ${IMAGE_TAG} run_tests"
    echo ""
    if [[ "${VERSION}" == @(sdbg|ssmp) ]]; then
      echo "*** A CP2K run using 2 OpenMP threads (default) can be launched with"
      echo "    podman run -it${DEVICE_FLAG} --rm ${IMAGE_TAG} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    else
      echo "*** An MPI-parallel CP2K run using 2 OpenMP threads for each of the 4 MPI ranks can be launched with"
      echo "    podman run -it${DEVICE_FLAG} --rm ${IMAGE_TAG} mpiexec -n 4 ${ENV_VAR_FLAG} OMP_NUM_THREADS=2 cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    fi
  else
    echo ""
    echo "*** A regression test run can be launched with"
    echo "    ${LAUNCH_SCRIPT} run_tests"
    echo ""
    if [[ "${VERSION}" == @(sdbg|ssmp) ]]; then
      echo "*** A CP2K run using 2 OpenMP threads (default) can be launched with"
      echo "    ${LAUNCH_SCRIPT} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
      echo ""
      echo "*** A CP2K run using only 4 OpenMP threads can be launched with"
      echo "    export OMP_NUM_THREADS=4; ${LAUNCH_SCRIPT} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    else
      echo "*** An MPI-parallel CP2K run using 2 OpenMP threads for each of the 4 MPI ranks can be launched with"
      echo "    export OMP_NUM_THREADS=2; ${LAUNCH_SCRIPT} mpiexec -n 4 cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    fi
  fi
fi

# Optionally, analyse code coverage and generate a coverage report
if [[ "${TEST_COVERAGE}" == "yes" ]]; then
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
fi

# Optionally, run CP2K benchmark as performance check
if [[ "${VERSION}" == "psmp" ]]; then
  if [[ "${RUN_BENCHMARK}" == "yes" ]]; then
    echo -e "\n*** Launching benchmark run using the script ${INSTALL_PREFIX}/bin/run_benchmarks\n"
    ${LAUNCH_SCRIPT} run_benchmarks
    EXIT_CODE=$?
    if ((EXIT_CODE != 0)); then
      echo "ERROR: The benchmark run failed with the error code ${EXIT_CODE}"
      exit "${EXIT_CODE}"
    fi
  else
    if [[ "${IN_CONTAINER}" != "yes" ]]; then
      echo ""
      echo "*** A benchmark run can be launched with"
      echo "    ${LAUNCH_SCRIPT} run_benchmarks"
      echo ""
    fi
  fi
fi
