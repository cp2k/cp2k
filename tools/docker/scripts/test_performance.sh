#!/bin/bash -e

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: test_performance.sh <ARCH>"
  exit 1
fi

ARCH=$1

function run_benchmark {
  set +e #disable error trapping
  OMP_THREADS=$1
  MPI_RANKS=$2
  INPUT=$3
  OUTPUT=$4
  echo -n "Running ${INPUT} with ${OMP_THREADS} threads and ${MPI_RANKS} ranks... "
  if OMP_NUM_THREADS="${OMP_THREADS}" mpiexec -np "${MPI_RANKS}" \
    "/opt/cp2k/exe/${ARCH}/cp2k.psmp" "${INPUT}" &> "${OUTPUT}"; then
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

echo -e '\n========== Compiling CP2K =========='
cd /opt/cp2k
echo -n "Compiling cp2k... "
if make -j ARCH="${ARCH}" VERSION="psmp" &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  cp make.out /workspace/artifacts/
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

# Check benchmark files for input errors.
echo -en "\nChecking benchmark inputs... "
if ! ./tools/regtesting/check_inputs.py "./exe/${ARCH}/cp2k.psmp" "./benchmarks/"; then
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

if [[ "${ARCH}" == "local" ]]; then
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
      "${LABEL} with 32 OpenMP Threads" "${LABEL}_timings_32omp" "${OUTPUT_OMP}" \
      "${LABEL} with 32 MPI Ranks" "${LABEL}_timings_32mpi" "${OUTPUT_MPI}"
    echo ""
  done

elif [[ "${ARCH}" == "local_cuda" ]]; then
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
      "${LABEL} with 6 CPU Cores and 1 GPU" "${LABEL}_timings_6cpu_1gpu" "${OUTPUT}"
    echo ""
  done

else
  echo "Unknown ARCH: ${ARCH}"
  exit 1
fi

TIME_END=$(date +%s)
DURATION=$(printf "%i" $(((TIME_END - TIME_START) / 60)))

echo -e "\nSummary: Performance test took ${DURATION} minutes."
echo -e "Status: OK\n"

#EOF
