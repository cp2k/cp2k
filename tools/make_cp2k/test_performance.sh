#!/usr/bin/env bash
#
# test_performance.sh: Create the CP2K benchmark script and run the benchmarks (step 4 of make_cp2k.sh)
#
# Usage: test_performance.sh CONFIG_FILE...
#
# The settings are read from the configuration files options.conf, deps.conf,
# and build.conf written by make_cp2k.sh and build_deps.sh.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"
activate_spack_env

# Create script to run the CP2K benchmarks for psmp builds
if [[ "${VERSION}" == "psmp" ]]; then
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    BENCHMARK_OUTPUT_DIR="/workspace/artifacts"
  else
    BENCHMARK_OUTPUT_DIR="${INSTALL_PREFIX}"/performance_tests
  fi
  cat << *** > "${INSTALL_PREFIX}"/bin/run_benchmarks
#!/bin/bash -e

BENCHMARK_PROFILE=\${1:-openmp}

MAX_PROCS=\$(nproc --all)
if ((MAX_PROCS < 32)); then
  echo -e "\nWARNING: 32 CPU cores are needed for a CP2K benchmark run but if seems only \${MAX_PROCS} are available"
fi

run_benchmark_input() {
  set +e
  local omp_threads="\$1"
  local mpi_ranks="\$2"
  local input="\$3"
  local output="\$4"
  printf "Running %s with %s threads and %s ranks... " \
    "\${input}" "\${omp_threads}" "\${mpi_ranks}"
  if OMP_NUM_THREADS="\${omp_threads}" mpiexec -n "\${mpi_ranks}" cp2k.psmp "\${input}" >"\${output}" 2>&1; then
    echo "done."
  else
    echo "failed."
    echo
    [[ -f "\${output}" ]] && tail -n 1000 "\${output}"
    echo
    echo "Summary: Running \${input} failed."
    echo "Status: FAILED"
    exit 0
  fi
  set -e
}

echo -e "\n========== Running Performance Test =========="
OUTPUT_DIR=${BENCHMARK_OUTPUT_DIR}
mkdir -p \${OUTPUT_DIR}
TIME_START=\$(date +%s)

BENCHMARKS=(
  "${CP2K_ROOT}/benchmarks/QS/H2O-64.inp"
  "${CP2K_ROOT}/benchmarks/QS/H2O-64_nonortho.inp"
  "${CP2K_ROOT}/benchmarks/QS_reference/w64PBE.inp"
  "${CP2K_ROOT}/benchmarks/QS_reference/w64GAPW-SVP.inp"
  "${CP2K_ROOT}/benchmarks/QS_reference/w64SCAN.inp"
  "${CP2K_ROOT}/benchmarks/QS_kp/ZnO.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/H2O-hyb.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/GW_PBE_4benzene.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/RI-HFX_H2O-32.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/RI-MP2_ammonia.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/diag_cu144_broy.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/bench_dftb.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/dbcsr.inp"
  "${CP2K_ROOT}/benchmarks/QMMM/MQAE/MQAE_single_node.inp"
)

if [[ "\${BENCHMARK_PROFILE}" == "openmp" ]]; then
  echo 'Plot: name="total_timings_32omp", title="Total Timings with 32 OpenMP Threads", ylabel="time [s]"'
  echo 'Plot: name="total_timings_32mpi", title="Total Timings with 32 MPI Ranks", ylabel="time [s]"'
  echo ""

  for INPUT in "\${BENCHMARKS[@]}"; do
    INPUT_BASENAME=\$(basename "\${INPUT}")
    LABEL=\${INPUT_BASENAME%.*}
    OUTPUT_MPI="\${OUTPUT_DIR}/\${LABEL}_32mpi.out"
    OUTPUT_OMP="\${OUTPUT_DIR}/\${LABEL}_32omp.out"
    cd "\$(dirname "\${INPUT}")"
    run_benchmark_input 1 32 "\${INPUT_BASENAME}" "\${OUTPUT_MPI}"
    run_benchmark_input 32 1 "\${INPUT_BASENAME}" "\${OUTPUT_OMP}"
    cd ..
    echo ""
    ${CP2K_ROOT}/tools/docker/scripts/plot_performance.py \
      "\${LABEL} with 32 OpenMP Threads" "\${LABEL}" "32omp" "\${OUTPUT_OMP}" \
      "\${LABEL} with 32 MPI Ranks" "\${LABEL}" "32mpi" "\${OUTPUT_MPI}"
    echo ""
  done

else
  echo "ERROR: Unknown benchmark profile \${BENCHMARK_PROFILE} specified"
  exit 1
fi

TIME_END=\$(date +%s)
DURATION=\$(printf "%i" \$(((TIME_END - TIME_START) / 60)))

echo -e "\nSummary: Performance test took \${DURATION} minutes."
echo -e "Status: OK\n"
***
  chmod 750 "${INSTALL_PREFIX}"/bin/run_benchmarks
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
