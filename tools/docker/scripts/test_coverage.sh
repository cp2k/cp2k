#!/bin/bash

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Get CP2K revision.
cd /opt/cp2k || exit 1
CP2K_REVISION=$(./tools/build_utils/get_revision_number ./src)

echo -e "\n========== Installing Dependencies =========="
lcov_version="2.3.1"
lcov_sha256="b3017679472d5fcca727254493d0eb44253c564c2c8384f86965ba9c90116704"

# LCOV dependencies
apt-get update -qq
apt-get install -qq --no-install-recommends \
  libperlio-gzip-perl \
  libcapture-tiny-perl \
  libdatetime-perl \
  libtimedate-perl \
  libjson-xs-perl
rm -rf /var/lib/apt/lists/*

cd /tmp || exit 1
wget -q "https://www.cp2k.org/static/downloads/lcov-${lcov_version}.tar.gz"
echo "${lcov_sha256} lcov-${lcov_version}.tar.gz" | sha256sum --check
tar -xzf "lcov-${lcov_version}.tar.gz"
cd "lcov-${lcov_version}" || exit 1
make install > make.log 2>&1
cd .. || exit 1
rm -rf "lcov-${lcov_version}.tar.gz" "lcov-${lcov_version}"

echo -e "\n========== Running Regtests =========="
cd /opt/cp2k || exit 1
./tests/do_regtest.py ./build/bin/ psmp

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
  --ignore-errors inconsistent \
  --capture \
  --output-file coverage.info > lcov.log

# print summary
lcov --summary coverage.info

# generate html report
genhtml \
  --ignore-errors range \
  --title "CP2K Regtests (${CP2K_REVISION})" \
  coverage.info > genhtml.log

# plot
LINE_COV=$(lcov --summary coverage.info | grep lines | awk '{print substr($2, 1, length($2)-1)}')
FUNC_COV=$(lcov --summary coverage.info | grep funct | awk '{print substr($2, 1, length($2)-1)}')
echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
echo "PlotPoint: name='lines', plot='cov', label='Lines', y=$LINE_COV, yerr=0"
echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=$FUNC_COV, yerr=0"

#EOF
