#!/bin/bash

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: test_coverage.sh <VERSION>"
  exit 1
fi

ARCH=local_coverage
VERSION=$1

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /opt/cp2k || exit 1
CP2K_REVISION=$(./tools/build_utils/get_revision_number ./src)
rm -rf "obj/${ARCH}/${VERSION}"/*.gcda # remove old gcov statistics

# Compile cp2k.
echo -en "\nCompiling cp2k... "
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Installing Dependencies =========="
lcov_version="1.15"
lcov_sha256="c1cda2fa33bec9aa2c2c73c87226cfe97de0831887176b45ee523c5e30f8053a"

# LCOV dependencies
apt-get update -qq
apt-get install -qq --no-install-recommends \
  libperlio-gzip-perl \
  libjson-perl
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
make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

# gcov gets stuck on some files...
# Maybe related: https://bugs.launchpad.net/gcc-arm-embedded/+bug/1694644
# As a workaround we'll simply remove the offending files for now.
rm -f "/opt/cp2k/obj/${ARCH}/${VERSION}"/exts/*/*.gcda
cd /tmp || exit 1
GCOV_TIMEOUT="10s"
for fn in "/opt/cp2k/obj/${ARCH}/${VERSION}"/*.gcda; do
  if ! timeout "${GCOV_TIMEOUT}" gcov "$fn" &> /dev/null; then
    echo "Skipping ${fn} because gcov took longer than ${GCOV_TIMEOUT}."
    rm "${fn}"
  fi
done

# collect coverage stats
mkdir -p /workspace/artifacts/coverage
cd /workspace/artifacts/coverage || exit 1
lcov --directory "/opt/cp2k/obj/${ARCH}/${VERSION}" \
  --exclude "/opt/cp2k-toolchain/*" \
  --capture \
  --output-file coverage.info > lcov.log
lcov --summary coverage.info

# generate html report
genhtml --title "CP2K Regtests (${CP2K_REVISION})" coverage.info > genhtml.log

# plot
LINE_COV=$(lcov --summary coverage.info | grep lines | awk '{print substr($2, 1, length($2)-1)}')
FUNC_COV=$(lcov --summary coverage.info | grep funct | awk '{print substr($2, 1, length($2)-1)}')
echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
echo "PlotPoint: name='lines', plot='cov', label='Lines', y=$LINE_COV, yerr=0"
echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=$FUNC_COV, yerr=0"

#EOF
