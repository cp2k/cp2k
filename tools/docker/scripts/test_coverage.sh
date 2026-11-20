#!/bin/bash

# author: Ole Schuett

ulimit -c 0 # Disable core dumps as they can take a very long time to write.

# Check available shared memory - needed for MPI inter-process communication.
SHM_AVAIL=$(df --output=avail -m /dev/shm | tail -1)
if ((SHM_AVAIL < 1024)); then
  echo "ERROR: Not enough shared memory. If you're running docker use --shm-size=1g."
  exit 1
fi

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Get CP2K revision.
cd /opt/cp2k || exit 1
CP2K_REVISION=$(./tools/build_utils/get_revision_number ./src)

echo -e "\n========== Running Regtests =========="
./tests/do_regtest.py --ompthreads=1 --keepalive ./build/bin/ psmp

echo -e "\n========== Analysing Coverage =========="
./tools/make_cp2k/test_coverage.sh \
  --build_dir /opt/cp2k/build/src \
  --output_dir /workspace/artifacts/coverage \
  --exclude "/opt/cp2k-toolchain/*" \
  --title "CP2K Regtests (${CP2K_REVISION})" \
  --gcov_timeout 10s

#EOF
