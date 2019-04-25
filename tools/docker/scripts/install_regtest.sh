#!/bin/bash

# author: Ole Schuett

if (( $# != 2 )) ; then
    echo "Usage: install_regtest.sh <ARCH> <VERSION>"
    exit 1
fi

ARCH=$1
VERSION=$2

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# pre-build cp2k
cd /workspace/cp2k
echo -n "Warming cache by trying to compile... "
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> /dev/null ; then
    echo "done."
else
    echo "failed."
fi


# run regtests which lack fixed reference value
# Disable LeakSanitizer during docker build as it requires ptrace capabilities.
export LSAN_OPTIONS="detect_leaks=0"
echo -n "Trying to run regtests which lack reference values... "
if make test ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="-restrictdir QS/regtest-almo-md -restrictdir QS/regtest-almo-1 -restrictdir SE/regtest-3-4 -restrictdir QS/regtest-ot-1-vib -restrictdir Fist/regtest-5-vib -restrictdir QS/regtest-optbas -restrictdir TMC/regtest_ana_post_proc" &> /dev/null ; then
   echo "done."
else
   echo "failed."
fi

# remove binaries to reduce image size
rm -rf lib exe "regtesting/${ARCH}/${VERSION}"/TEST-*

#EOF
