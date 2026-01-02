#!/bin/bash -e

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: test_performance.sh <PROFILE>"
  exit 1
fi

PROFILE=$1

function run_benchmark {
  set +e #disable error trapping
  OMP_THREADS=$1
  MPI_RANKS=$2
  INPUT=$3
  OUTPUT=$4
  echo -n "Running ${INPUT} with ${OMP_THREADS} threads and ${MPI_RANKS} ranks... "
  if OMP_NUM_THREADS="${OMP_THREADS}" mpiexec -np "${MPI_RANKS}" \
    "/opt/cp2k/build/bin/cp2k.psmp" "${INPUT}" &> "${OUTPUT}"; then
    echo "done."
  else
    echo -e "failed.\n\n"
    tail -n 100 "${OUTPUT}"
    echo -e "\nSummary: Running ${INPUT} failed."
    echo -e "Status: FAILED\n"
    exit 0
  fi
  set -e #re-enable error trapping
}

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e '\n============== CP2K Binary Flags ============='
./build/bin/cp2k.psmp --version | grep cp2kflags

echo -e '\n========== Checking Benchmark Inputs ========='
if ! ./tools/regtesting/check_inputs.py "./build/bin/cp2k.psmp" "./benchmarks/"; then
  echo -e "\nSummary: Some benchmark inputs could not be parsed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e '\n========== Running Performance Test =========='
mkdir -p /workspace/artifacts
cd ./benchmarks
TIME_START=$(date +%s)

BENCHMARKS=(
  "QS/H2O-64.inp"
  "QS/H2O-64_nonortho.inp"
  "QS_single_node/H2O-hyb.inp"
  "QS_single_node/GW_PBE_4benzene.inp"
  "QS_single_node/RI-HFX_H2O-32.inp"
  "QS_single_node/RI-MP2_ammonia.inp"
  "QS_single_node/diag_cu144_broy.inp"
  "QS_single_node/bench_dftb.inp"
  "QS_single_node/dbcsr.inp"
  "QMMM_MQAE/MQAE_single_node.inp"
)

if [[ "${PROFILE}" == "toolchain" ]]; then
  echo 'Plot: name="total_timings_32omp", title="Total Timings with 32 OpenMP Threads", ylabel="time [s]"'
  echo 'Plot: name="total_timings_32mpi", title="Total Timings with 32 MPI Ranks", ylabel="time [s]"'
  echo ''

  for INPUT in "${BENCHMARKS[@]}"; do
    INPUT_BASENAME=$(basename "${INPUT}")
    LABEL=${INPUT_BASENAME%.*}
    OUTPUT_MPI="/workspace/artifacts/${LABEL}_32mpi.out"
    OUTPUT_OMP="/workspace/artifacts/${LABEL}_32omp.out"
    cd "$(dirname "${INPUT}")"
    run_benchmark 1 32 "${INPUT_BASENAME}" "${OUTPUT_MPI}"
    run_benchmark 32 1 "${INPUT_BASENAME}" "${OUTPUT_OMP}"
    cd ..
    echo ""
    /opt/cp2k/plot_performance.py \
      "${LABEL} with 32 OpenMP Threads" "${LABEL}" "32omp" "${OUTPUT_OMP}" \
      "${LABEL} with 32 MPI Ranks" "${LABEL}" "32mpi" "${OUTPUT_MPI}"
    echo ""
  done

elif [[ "${PROFILE}" == "toolchain_cuda_"* ]]; then
  echo 'Plot: name="total_timings_6cpu_1gpu", title="Total Timings with 6 CPU Cores and 1 GPU", ylabel="time [s]"'
  echo ''

  for INPUT in "${BENCHMARKS[@]}"; do
    if [[ "$INPUT" == "QS_single_node/H2O-hyb.inp" ]]; then
      continue # Has no gpu acceleration, yet.
    fi
    INPUT_BASENAME=$(basename "${INPUT}")
    LABEL=${INPUT_BASENAME%.*}
    OUTPUT="/workspace/artifacts/${LABEL}_6cpu_1gpu.out"
    cd "$(dirname "${INPUT}")"
    export CUDA_VISIBLE_DEVICES=0
    run_benchmark 3 2 "${INPUT_BASENAME}" "${OUTPUT}"
    cd ..
    echo ""
    /opt/cp2k/plot_performance.py \
      "${LABEL} with 6 CPU Cores and 1 GPU" "${LABEL}" "6cpu_1gpu" "${OUTPUT}"
    echo ""
  done

else
  echo "Unknown PROFILE: ${PROFILE}"
  exit 1
fi

TIME_END=$(date +%s)
DURATION=$(printf "%i" $(((TIME_END - TIME_START) / 60)))

echo -e "\nSummary: Performance test took ${DURATION} minutes."
echo -e "Status: OK\n"

#EOF
