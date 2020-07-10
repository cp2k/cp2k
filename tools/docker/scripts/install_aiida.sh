#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update -qq
export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true
apt-get install -qq --no-install-recommends \
    python3-setuptools    \
    python3-wheel         \
    python3-pip           \
    python3-dev           \
    postgresql            \
    rabbitmq-server       \
    sudo                  \
    ssh
rm -rf /var/lib/apt/lists/*

# install dependencies of aiida-cp2k
cd /opt/
git clone --quiet https://github.com/aiidateam/aiida-cp2k.git
pip3 install --quiet './aiida-cp2k/[test]'
pip3 uninstall --quiet --yes aiida-cp2k

# create ubuntu user with sudo powers
adduser --disabled-password --gecos "" ubuntu
echo "ubuntu ALL=(ALL) NOPASSWD: ALL" >>  /etc/sudoers

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# link mpi executables into path
MPI_INSTALL_DIR=$(dirname "$(command -v mpirun)")
for i in "${MPI_INSTALL_DIR}"/* ; do ln -sf "$i" /usr/bin/; done

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
