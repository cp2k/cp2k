#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true
apt-get install -y --no-install-recommends \
    build-essential       \
    python-setuptools     \
    python-wheel          \
    python-pip            \
    python-dev            \
    postgresql            \
    sudo                  \
    ssh
rm -rf /var/lib/apt/lists/*

# install python packages
pip install flake8 aiida ase

# create ubuntu user with sudo powers
adduser --disabled-password --gecos "" ubuntu
echo "ubuntu ALL=(ALL) NOPASSWD: ALL" >>  /etc/sudoers

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
