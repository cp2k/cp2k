#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true
apt-get install -qq --no-install-recommends \
    build-essential       \
    python-setuptools     \
    python-wheel          \
    python-pip            \
    python-dev            \
    postgresql            \
    rabbitmq-server       \
    sudo                  \
    ssh
rm -rf /var/lib/apt/lists/*

# install python packages
pip install --quiet flake8 aiida ase

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
echo -n "Warming cache by trying to compile... "
if make -j VERSION=pdbg &> /dev/null ; then
   echo "done."
else
   echo "failed."
fi
rm -rf lib exe

#EOF
