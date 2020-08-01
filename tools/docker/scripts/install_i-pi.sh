#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
apt-get install -qq --no-install-recommends \
    python3                                 \
    python3-pip                             \
    python3-wheel                           \
    python3-setuptools
rm -rf /var/lib/apt/lists/*

# install python packages
pip3 install --quiet numpy

# clone i-pi repository
git clone --quiet --depth=1 --single-branch -b master https://github.com/i-pi/i-pi.git /opt/i-pi

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# pre-build cp2k
cd /workspace/cp2k
echo -n "Warming cache by trying to compile... "
if make -j VERSION=pdbg &> /dev/null ; then
   echo "done."
else
   echo "failed."
fi
rm -rf lib exe

#EOF
