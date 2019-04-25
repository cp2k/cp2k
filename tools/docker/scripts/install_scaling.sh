#!/bin/bash -e

# author: Ole Schuett

# install numpy
apt-get update -qq
apt-get install -qq --no-install-recommends python-numpy
rm -rf /var/lib/apt/lists/*

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# pre-build cp2k
cd /workspace/cp2k
echo -n "Warming cache by trying to compile... "
if make -j VERSION="popt" cp2k &> /dev/null && \
   make -j VERSION="psmp" cp2k &> /dev/null && \
   make -j VERSION="ssmp" cp2k &> /dev/null ; then
   echo "done."
else
   echo "failed."
fi


rm -rf lib exe

#EOF
