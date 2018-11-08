#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends  \
    default-jre-headless                    \
    libsaxonhe-java
rm -rf /var/lib/apt/lists/*

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# pre-build cp2k
cd /workspace/cp2k
make -j VERSION=sopt
rm -rf lib exe

#EOF
