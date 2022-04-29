#!/bin/bash

# author: Ole Schuett

if (($# != 2)); then
  echo "Usage: test_regtest.sh <ARCH> <VERSION>"
  exit 1
fi

ARCH=$1
VERSION=$2

ulimit -c 0 # Disable core dumps as they can take a very long time to write.

# Extend stack size - needed when using Intel compilers.
ulimit -s unlimited
export OMP_STACKSIZE=64m

# Check available shared memory - needed for MPI inter-process communication.
SHM_AVAIL=$(df --output=avail -m /dev/shm | tail -1)
if ((SHM_AVAIL < 1024)); then
  echo "ERROR: Not enough shared memory. If you're running docker use --shm-size=1g."
  exit 1
fi

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Make OpenMPI happy.
if command -v ompi_info &> /dev/null; then
  TESTOPTS="--mpiexec='mpiexec --bind-to none --allow-run-as-root' ${TESTOPTS}"
  export OMPI_MCA_plm_rsh_agent=/bin/false
fi

# Switch to stable DBCSR version if requested.
if [ -n "${USE_STABLE_DBCSR}" ]; then
  echo "Switching to stable DBCSR version..."
  if ! git -C cp2k/exts/dbcsr checkout v2.1.0-rc16; then
    echo -e "\nSummary: Could not checkout stable DBCSR version."
    echo -e "Status: FAILED\n"
    exit 0
  fi
  ln -fs python3 /usr/bin/python # DBCSR v2.1.0-rc16 needs the python binary.
fi

# Compile cp2k.
echo -en "\nCompiling cp2k... "
cd /opt/cp2k || exit 1
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

# Improve code coverage on COSMA.
export COSMA_DIM_THRESHOLD=0

# Run regtests.
echo -e "\n========== Running Regtests =========="
make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

exit 0 # Prevent CI from overwriting do_regtest's summary message.

#EOF
