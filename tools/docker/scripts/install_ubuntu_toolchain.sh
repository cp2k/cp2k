#!/bin/bash -e

# author: Ole Schuett

if (( $# != 1 )) ; then
    echo "Usage: install_ubuntu_toolchain.sh <GCC_VERSION>"
    exit 1
fi

GCC_VERSION=$1

# install Ubuntu packages
apt-get update
apt-get install -y --no-install-recommends  \
    gcc-${GCC_VERSION}                      \
    g++-${GCC_VERSION}                      \
    gfortran-${GCC_VERSION}                 \
    fftw3-dev                               \
    libopenblas-dev                         \
    liblapack-dev                           \
    libint-dev                              \
    libgsl-dev                              \
    libhdf5-dev
rm -rf /var/lib/apt/lists/*

# create links
ln -sf gcc-${GCC_VERSION}      /usr/bin/gcc
ln -sf g++-${GCC_VERSION}      /usr/bin/g++
ln -sf gfortran-${GCC_VERSION} /usr/bin/gfortran

# json-fortran does not compile with gcc 4.8.5.
if [[ "$GCC_VERSION" == "4.8" ]] ; then
  EXTRA_TOOLCHAIN_OPTIONS="--with-sirius=no --with-json-fortran=no"
fi

# build toolchain relying mostly on ubuntu packages
cp -r /workspace/cp2k/tools/toolchain /opt/cp2k-toolchain/
cd /opt/cp2k-toolchain/
./install_cp2k_toolchain.sh  \
    --mpi-mode=no            \
    --with-gcc=system        \
    --with-binutils=system   \
    --with-make=system       \
    --with-cmake=system      \
    --with-fftw=system       \
    --with-openblas=system   \
    --with-reflapack=system  \
    --with-libint=system     \
    --with-gsl=system        \
    --with-hdf5=system       \
    --with-libxc=install     \
    --with-libxsmm=install   \
    ${EXTRA_TOOLCHAIN_OPTIONS}
rm -rf ./build

#EOF
