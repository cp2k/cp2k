#!/bin/bash -e

# author: Ole Schuett

if (( $# != 1 )) ; then
    echo "Usage: test_coverage.sh <VERSION>"
    exit 1
fi

ARCH=local_coverage
VERSION=$1

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Running Regtests =========="
cd /workspace/cp2k
CP2K_REVISION=$(./tools/build_utils/get_revision_number)
rm -rf "obj/${ARCH}/${VERSION}"/*.gcda   # remove old gcov statistics

make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

# gcov gets stuck on some files...
# Maybe related: https://bugs.launchpad.net/gcc-arm-embedded/+bug/1694644
# As a workaround we'll simply remove the offending files for now.
rm -v "/workspace/cp2k/obj/${ARCH}/${VERSION}"/almo_scf_types.gcda
rm -v "/workspace/cp2k/obj/${ARCH}/${VERSION}"/dbcsr_tensor_types.gcda
rm -v "/workspace/cp2k/obj/${ARCH}/${VERSION}"/mp2_types.gcda
#cd /tmp
#for i in /workspace/cp2k/obj/${ARCH}/${VERSION}/*.gcda; do
#    timeout 30 gcov $i >/dev/null 2>&1 || rm -v $i
#done

# collect coverage stats
mkdir -p /workspace/artifacts/coverage
cd /workspace/artifacts/coverage
lcov --directory "/workspace/cp2k/obj/${ARCH}/${VERSION}" --capture --output-file coverage.info > lcov.log
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
