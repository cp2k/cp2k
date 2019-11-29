#!/bin/bash

# author: Ole Schuett

if (( $# != 1 )) ; then
    echo "Usage: test_benchmarks.sh <ARCH>"
    exit 1
fi

ARCH=$1

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

cd /workspace/cp2k
echo -n "Compiling cp2k... "
if make -j ARCH="${ARCH}" VERSION="psmp" &> /dev/null ; then
   echo "done."
else
   echo -e "failed.\n\n"
   echo "Summary: Compilation failed."
   echo "Status: FAILED"
   exit
fi

echo -e "\n========== Running Benchmarks =========="
ARTIFACTS="/workspace/artifacts/benchmarks/"
mkdir -p ${ARTIFACTS}

NPROC=$(nproc) # Call nproc before OMP_NUM_THREADS is exported to get actual number of procs.
export OMP_NUM_THREADS=8
NUM_MPI_RANKS=$(( NPROC / OMP_NUM_THREADS ))
echo "Found ${NPROC} processors, using ${NUM_MPI_RANKS} MPI ranks with ${OMP_NUM_THREADS} OpenMP threads each."
CP2K="mpiexec -np ${NUM_MPI_RANKS} /workspace/cp2k/exe/${ARCH}/cp2k.psmp"

echo "Running benchmarks/QS/H2O-64.inp ..."
T1=$(date +%s)
cd /workspace/cp2k/benchmarks/QS
${CP2K} H2O-64.inp >& "${ARTIFACTS}H2O-64.out"
T2=$(date +%s)
H2O_OT_TIMING=$((1 + T2 - T1))

echo "Running benchmarks/QS_DM_LS/H2O-dft-ls.NREP2.inp ..."
T1=$(date +%s)
cd /workspace/cp2k/benchmarks/QS_DM_LS/
${CP2K} H2O-dft-ls.NREP2.inp >& "${ARTIFACTS}H2O-dft-ls.NREP2.out"
T2=$(date +%s)
H2O_LS_TIMING=$((1 + T2 - T1))

echo ""
echo "Plot: name='timings', title='Benchmark Timings', ylabel='time [s]'"
echo "PlotPoint: plot='timings', name='h2o-ot-64',    label='H2O-64 MD',    yerr=0, y=${H2O_OT_TIMING}"
echo "PlotPoint: plot='timings', name='h2o-ls-nrep2', label='H2O LS NREP2', yerr=0, y=${H2O_LS_TIMING}"
echo ""
echo "Summary: H2O-OT: ${H2O_OT_TIMING}s, H2O-LS: ${H2O_LS_TIMING}s"
echo "Status: OK"
echo ""

#EOF
