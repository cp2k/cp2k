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

   mkdir -p /opt/cp2k_test_artifacts/coverage
   cd /opt/cp2k_test_artifacts/coverage
   lcov --directory /opt/cp2k-master/cp2k/obj/${ARCH}/${VERSION} --capture --output-file coverage.info > lcov.log
   genhtml coverage.info > genhtml.log

   # plot
   DA=$(grep "^DA:" coverage.info  | wc -l)
   DA_ZERO=$(grep "^DA:.*,0" coverage.info  | wc -l)
   DA_PERCENT=$(python -c "print(100.0 * ($DA - $DA_ZERO) / $DA)")
   FNDA=$(grep "^FNDA:" coverage.info  | wc -l)
   FNDA_ZERO=$(grep "^FNDA:0" coverage.info  | wc -l)
   FNDA_PERCENT=$(python -c "print(100.0 * ($FNDA - $FNDA_ZERO) / $FNDA)")
   echo 'Plot: name="cov", title="Test Coverage", ylabel="Coverage %"'
   echo "PlotPoint: name='lines', plot='cov', label='Lines', y=$DA_PERCENT, yerr=0"
   echo "PlotPoint: name='funcs', plot='cov', label='Functions', y=$FNDA_PERCENT, yerr=0"
fi
#EOF
