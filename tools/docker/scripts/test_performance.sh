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
cd ./benchmarks/QS

for INPUT in "H2O-64.inp" "H2O-64_nonortho.inp"; do
  LABEL="${INPUT%.*}"
  OUTPUT_MPI="/workspace/artifacts/${LABEL}_32mpi.out"
  OUTPUT_OMP="/workspace/artifacts/${LABEL}_32omp.out"
  run_benchmark 1 32 "${INPUT}" "${OUTPUT_MPI}"
  run_benchmark 32 1 "${INPUT}" "${OUTPUT_OMP}"
  echo ""
  /workspace/plot_performance.py "${LABEL}" "${LABEL}" "${OUTPUT_OMP}" "${OUTPUT_MPI}"
  echo ""
done

echo -e "\nSummary: Performance test works fine."
echo -e "Status: OK\n"

#EOF
