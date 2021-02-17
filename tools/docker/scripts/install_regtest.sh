#!/bin/bash

# author: Ole Schuett

if (($# != 2)); then
  echo "Usage: install_regtest.sh <ARCH> <VERSION>"
  exit 1
fi

ARCH=$1
VERSION=$2

# setup arch files
cd /workspace/cp2k/arch || exit 1
ln -vs /opt/cp2k-toolchain/install/arch/local* .

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Make OpenMPI happy.
if command -v ompi_info &> /dev/null; then
  TESTOPTS="-mpiexec 'mpiexec --bind-to none --allow-run-as-root' ${TESTOPTS}"
  export OMPI_MCA_plm_rsh_agent=/bin/false
fi

# pre-build cp2k
cd /workspace/cp2k || exit 1
echo -n "Warming cache by trying to compile... "
if make -j ARCH="${ARCH}" VERSION="${VERSION}" &> /dev/null; then
  echo "done."
else
  echo "failed."
fi

# remove binaries to reduce image size
rm -rf lib exe "regtesting/${ARCH}/${VERSION}"/TEST-*

#EOF
