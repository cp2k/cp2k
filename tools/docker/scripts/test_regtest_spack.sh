#!/bin/bash

ulimit -c 0 # Disable core dumps as they can take a very long time to write.

# Check available shared memory - needed for MPI inter-process communication.
SHM_AVAIL=$(df --output=avail -m /dev/shm | tail -1)
if ((SHM_AVAIL < 1024)); then
  echo "ERROR: Not enough shared memory. If you're running docker use --shm-size=1g."
  exit 1
fi

# Improve code coverage on COSMA.
export COSMA_DIM_THRESHOLD=0

ulimit -s unlimited
export OMP_STACKSIZE=64m

# Run regtests.
echo -e "\n========== Running Regtests =========="
set -x
./tests/do_regtest.py local psmp ${TESTOPTS:+""}

exit 0 # Prevent CI from overwriting do_regtest's summary message.
