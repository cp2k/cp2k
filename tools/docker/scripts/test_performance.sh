#!/bin/bash -e

# author: Ole Schuett

function run_benchmark {
  set +e #disable error trapping
  OMP_THREADS=$1
  MPI_RANKS=$2
  INPUT=$3
  OUTPUT=$4
  echo -n "Running ${INPUT} with ${OMP_THREADS} threads and ${MPI_RANKS} ranks... "
  if OMP_NUM_THREADS="${OMP_THREADS}" mpiexec -np "${MPI_RANKS}" \
    /workspace/cp2k/exe/local/cp2k.psmp "${INPUT}" &> "${OUTPUT}"; then
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
cd /workspace/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION="psmp" &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e '\n========== Running Performance Test =========='
mkdir -p /workspace/artifacts
cd ./benchmarks

BENCHMARKS=(
  "QS/H2O-64.inp"
  "QS/H2O-64_nonortho.inp"
  "QS_single_node/H2O-hyb.inp"
  "QS_single_node/bench_dftb.inp"
  "QS_single_node/dbcsr.inp"
)

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
  /workspace/plot_performance.py 32 "${LABEL}" "${LABEL}" "${OUTPUT_OMP}" "${OUTPUT_MPI}"
  echo ""
done

echo -e "\nSummary: Performance test works fine."
echo -e "Status: OK\n"

#EOF
