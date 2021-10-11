#!/bin/bash -e

# author: Ole Schuett

# clone gromacs reprository
git clone --quiet --depth=1 --single-branch -b master https://gitlab.com/gromacs/gromacs.git /opt/gromacs

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# setup arch files
cd /workspace/cp2k/arch
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# pre-build libcp2k
cd /workspace/cp2k
echo -n "Warming cache by trying to compile libcp2k... "
if make -j VERSION=pdbg libcp2k &> /dev/null; then
  echo "done."
else
  echo "failed."
fi
rm -rf lib exe

#EOF
