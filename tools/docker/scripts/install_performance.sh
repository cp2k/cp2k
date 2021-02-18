#!/bin/bash -e

# author: Ole Schuett

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# pre-build cp2k
cd /workspace/cp2k
echo -n "Warming cache by trying to compile cp2k... "
if make -j VERSION="psmp" &> /dev/null; then
  echo 'done.'
else
  echo 'failed.'
fi

rm -rf lib exe

#EOF
