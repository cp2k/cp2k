#!/bin/bash

# author: Ole Schuett
set -e

echo -e "\n========== Copying Changed Files =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/cp2k/obj/    \
      --exclude=/cp2k/lib/    \
      --exclude=/cp2k/exe/    \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

echo -e "\n========== Updating Toolchain =========="
cd /opt/cp2k-master/cp2k/tools/toolchain
./install_cp2k_toolchain.sh \
    --mpi-mode=no \
    --with-gcc=system \
    --with-binutils=system \
    --with-make=system \
    --with-fftw=system \
    --with-openblas=system \
    --with-reflapack=system \
    --with-libint=system \
    --with-libxc=install \
    --with-libxsmm=install

echo -e "\n========== Running Regtests =========="
source /opt/cp2k-master/cp2k/tools/toolchain/install/setup
cd /opt/cp2k-master/cp2k/makefiles
make VERSION=ssmp test TESTOPTS="${TESTOPTS}"

#EOF
