#!/bin/bash

# author: Ole Schuett

if (($# != 2)); then
  echo "Usage: test_regtest_cmake.sh <PROFILE> <VERSION>"
  exit 1
fi

PROFILE=$1
VERSION=$2

ulimit -c 0 # Disable core dumps as they can take a very long time to write.

# Check available shared memory - needed for MPI inter-process communication.
SHM_AVAIL=$(df --output=avail -m /dev/shm | tail -1)
if ((SHM_AVAIL < 1024)); then
  echo "ERROR: Not enough shared memory. If you're running docker use --shm-size=1g."
  exit 1
fi

# Compile CP2K.
./build_cp2k_cmake.sh "${PROFILE}" "${VERSION}" || exit 0

# Fake installation of data files.
mkdir -p ./share/cp2k
ln -s ../../data ./share/cp2k/data

# Increase stack size.
ulimit -s unlimited
export OMP_STACKSIZE=64m

# Improve code coverage on COSMA.
export COSMA_DIM_THRESHOLD=0

# Bind pika threads to first two cores. This is a hack. Do not use for production!
export PIKA_PROCESS_MASK="0x3"

# Load Spack or Toolchain environment.
if [[ "${PROFILE}" == "spack" ]]; then
  eval "$(spack env activate myenv --sh)"
elif [[ "${PROFILE}" == "toolchain" ]]; then
  # shellcheck disable=SC1091
  source /opt/cp2k-toolchain/install/setup
fi

# Run regtests.
echo -e "\n========== Running Regtests =========="
set -x
# shellcheck disable=SC2086
./tests/do_regtest.py local ${VERSION} ${TESTOPTS}

exit 0 # Prevent CI from overwriting do_regtest's summary message.

#EOF
