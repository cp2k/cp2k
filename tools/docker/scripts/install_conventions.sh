#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# pre-build cp2k
cd /workspace/cp2k/tools/conventions
./test_conventions.sh
rm -rf ../lib/ ../exe/

#EOF
