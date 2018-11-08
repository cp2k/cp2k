#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends \
    python                                 \
    python-pip                             \
    python-wheel                           \
    python-setuptools
rm -rf /var/lib/apt/lists/*

# install python packages
pip install numpy

# clone i-pi repository
git clone --depth=1 --single-branch -b master https://github.com/i-pi/i-pi.git /opt/i-pi

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# pre-build cp2k
cd /workspace/cp2k
make -j VERSION=pdbg
rm -rf lib exe

#EOF
