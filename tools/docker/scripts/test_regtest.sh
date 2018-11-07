#!/bin/bash -e

# author: Ole Schuett

ARCH=$1
VERSION=$2

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Running Regtests =========="
cd /workspace/cp2k
make ARCH="${ARCH}" VERSION="${VERSION}" TESTOPTS="${TESTOPTS}" test

#EOF
