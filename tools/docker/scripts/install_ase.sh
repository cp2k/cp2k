#!/bin/bash -e

# author: Ole Schuett

# install python packages
apt-get update
apt-get install -y --no-install-recommends \
    python3                                                      \
    python3-dev                                                  \
    python3-pip                                                  \
    python3-wheel                                                \
    python3-setuptools                                           \
    build-essential
rm -rf /var/lib/apt/lists/*

# install python packages
pip3 install numpy scipy matplotlib flask

# clone ase reprository
git clone --depth=1 --single-branch -b master https://gitlab.com/ase/ase.git /opt/ase

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
