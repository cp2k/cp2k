#!/bin/bash -e

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

ln -s /opt/cp2k/tools/build_utils/fypp /bin/fypp

mkdir build
cd build

if ! cmake \
  -Werror=dev \
  -DCP2K_USE_VORI=ON \
  -DCP2K_USE_COSMA=NO \
  -DCP2K_BLAS_VENDOR=OpenBLAS \
  -DCP2K_USE_SPGLIB=ON \
  -DCP2K_USE_LIBINT2=ON \
  -DCP2K_USE_LIBXC=ON \
  -DCP2K_USE_LIBTORCH=OFF \
  -DCP2K_BUILD_DBCSR=ON \
  -DCP2K_ENABLE_REGTESTS=ON \
  .. |& tee ./cmake.log; then
  echo -e "\nSummary: CMake failed."
  echo -e "Status: FAILED\n"
fi

if grep -A5 'CMake Warning' ./cmake.log; then
  echo -e "\nSummary: Found CMake warnings."
  echo -e "Status: FAILED\n"
fi

echo -e '\n\nCompiling cp2k...'
if make -j; then
  echo -e "\nSummary: Compilation works fine."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
fi

#EOF
