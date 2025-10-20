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

echo -e "\n========== Installing lcov =========="
apt-get update -qq
apt-get install -qq --no-install-recommends lcov libjson-xs-perl
rm -rf /var/lib/apt/lists/*

echo -e "\n========== Running Regtests =========="
cd /opt/cp2k || exit 1
./tests/do_regtest.py --ompthreads=1 ./build/bin/ psmp

# gcov gets stuck on some files...
# Maybe related: https://bugs.launchpad.net/gcc-arm-embedded/+bug/1694644
# As a workaround we'll simply remove the offending files for now.
cd /tmp || exit 1
GCOV_TIMEOUT="10s"
find /opt/cp2k/build/src/ -name "*.gcda" -print0 | while read -r -d $'\0' fn; do
  if ! timeout "${GCOV_TIMEOUT}" gcov "$fn" &> /dev/null; then
    echo "Skipping ${fn} because gcov took longer than ${GCOV_TIMEOUT}."
    rm "${fn}"
  fi
done

mkdir -p /workspace/artifacts/coverage
cd /workspace/artifacts/coverage || exit 1

# collect coverage stats
lcov --directory "/opt/cp2k/build/src" \
  --exclude "/opt/cp2k-toolchain/*" \
  --exclude "/usr/*" \
  --capture \
  --keep-going \
  --output-file coverage.info &> lcov.log

# print summary
lcov --summary coverage.info

# generate html report
genhtml \
  --keep-going \
  --title "CP2K Regtests (${CP2K_REVISION})" \
  coverage.info &> genhtml.log

# plot
LINE_COV=$(lcov --summary coverage.info | grep lines | awk '{print substr($2, 1, length($2)-1)}')
FUNC_COV=$(lcov --summary coverage.info | grep funct | awk '{print substr($2, 1, length($2)-1)}')
echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
echo "PlotPoint: name='lines', plot='cov', label='Lines', y=$LINE_COV, yerr=0"
echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=$FUNC_COV, yerr=0"

#EOF
