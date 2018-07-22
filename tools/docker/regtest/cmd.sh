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

#EOF
