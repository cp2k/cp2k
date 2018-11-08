#!/bin/bash -e

# author: Ole Schuett

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends  \
    build-essential                         \
    gfortran                                \
    fftw3-dev                               \
    libopenblas-dev                         \
    liblapack-dev                           \
    libint-dev
rm -rf /var/lib/apt/lists/*

# build toolchain relying mostly on ubuntu packages
cp -r /workspace/cp2k/tools/toolchain /opt/cp2k-toolchain/
cd /opt/cp2k-toolchain/
./install_cp2k_toolchain.sh  \
    --mpi-mode=no            \
    --with-gcc=system        \
    --with-binutils=system   \
    --with-make=system       \
    --with-fftw=system       \
    --with-openblas=system   \
    --with-reflapack=system  \
    --with-libint=system     \
    --with-libxc=install     \
    --with-libxsmm=install
rm -rf ./build

#EOF
