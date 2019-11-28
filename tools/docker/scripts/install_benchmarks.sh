#!/bin/bash -e

# author: Ole Schuett

if (( $# != 1 )) ; then
    echo "Usage: install_benchmarks.sh <ARCH>"
    exit 1
fi

ARCH=$1

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# pre-build cp2k
cd /workspace/cp2k
echo -n "Warming cache by trying to compile cp2k... "
if make -j ARCH="${ARCH}" VERSION="psmp" &> /dev/null ; then
    echo "done."
else
    echo "failed."
fi

# remove binaries to reduce image size
rm -rf lib exe

#EOF
