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

# Extend stack size only for Intel compilers.
if "./build/bin/cp2k.${VERSION}" --version | grep -q "compiler: Intel"; then
  ulimit -s unlimited # breaks address sanitizer
  export OMP_STACKSIZE=64m
fi

# Improve code coverage on COSMA.
export COSMA_DIM_THRESHOLD=0

# Make OpenMPI happy.
export OMPI_MCA_plm_rsh_agent=/bin/false
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Load Spack or Toolchain environment.
if [[ "${PROFILE}" =~ ^spack ]]; then
  eval "$(spack env activate myenv --sh)"
elif [[ "${PROFILE}" =~ ^toolchain ]]; then
  # shellcheck disable=SC1091
  source /opt/cp2k-toolchain/install/setup
fi

# Run regtests.
echo -e "\n========== Running Regtests =========="
set -x
# shellcheck disable=SC2086
./tests/do_regtest.py ./build/bin/ ${VERSION} ${TESTOPTS}

exit 0 # Prevent CI from overwriting do_regtest's summary message.

#EOF
