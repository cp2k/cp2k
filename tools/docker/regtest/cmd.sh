#!/bin/bash

# author: Ole Schuett
set -e

echo -e "\n========== Copying Changed Files =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/cp2k/obj/    \
      --exclude=/cp2k/lib/    \
      --exclude=/cp2k/exe/    \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

rsync --exclude="*~"                              \
      --exclude=".*/"                             \
      --exclude="*.pyc"                           \
      --ignore-times                              \
      --update                                    \
      --verbose                                   \
      --recursive                                 \
      --checksum                                  \
      /opt/cp2k-local/cp2k/tools/toolchain/  /opt/cp2k-toolchain/

echo -e "\n========== Updating Toolchain =========="
cd /opt/cp2k-toolchain/
./install_cp2k_toolchain.sh --install-all --with-make=no

echo -e "\n========== Running Regtests =========="
source /opt/cp2k-toolchain/install/setup
cd /opt/cp2k-master/cp2k/makefiles
rm -rf ../obj/${ARCH}/${VERSION}/*.gcda   # remove old gcov statistics

if [[ "$TESTNAME" != "farming" ]]; then
   make ARCH=${ARCH} VERSION=${VERSION} test TESTOPTS="${TESTOPTS}"
else
   make ARCH=${ARCH} VERSION=${VERSION} test TESTOPTS="-farming -skipunittest -skipdir TMC/regtest_ana_on_the_fly -skipdir TMC/regtest_ana_post_proc -skipdir TMC/regtest ${TESTOPTS}"
fi

if [[ "$TESTNAME" == coverage-* ]]; then
   # gcov gets stuck on some files...
   # Maybe related: https://bugs.launchpad.net/gcc-arm-embedded/+bug/1694644
   # As a workaround we'll simply remove the offending files for now.
   rm -v /opt/cp2k-master/cp2k/obj/${ARCH}/${VERSION}/almo_scf_types.gcda
   rm -v /opt/cp2k-master/cp2k/obj/${ARCH}/${VERSION}/dbcsr_tensor_types.gcda
   rm -v /opt/cp2k-master/cp2k/obj/${ARCH}/${VERSION}/mp2_types.gcda
   #cd /tmp
   #for i in /opt/cp2k-master/cp2k/obj/${ARCH}/${VERSION}/*.gcda; do
   #    timeout 30 gcov $i >/dev/null 2>&1 || rm -v $i
   #done

   # collect coverage stats
   mkdir -p /opt/cp2k_test_artifacts/coverage
   cd /opt/cp2k_test_artifacts/coverage
   lcov --directory /opt/cp2k-master/cp2k/obj/${ARCH}/${VERSION} --capture --output-file coverage.info > lcov.log
   lcov --summary coverage.info

   # generate html report
   genhtml --title "CP2K Regtests (${CP2K_REVISION})" coverage.info > genhtml.log

   # plot
   LINE_COV=$(lcov --summary coverage.info | grep lines | awk '{print substr($2, 1, length($2)-1)}')
   FUNC_COV=$(lcov --summary coverage.info | grep funct | awk '{print substr($2, 1, length($2)-1)}')
   echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
   echo "PlotPoint: name='lines', plot='cov', label='Lines', y=$LINE_COV, yerr=0"
   echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=$FUNC_COV, yerr=0"
fi
#EOF
