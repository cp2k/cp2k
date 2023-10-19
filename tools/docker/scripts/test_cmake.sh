#!/bin/bash -e

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

ln -s /opt/cp2k/tools/build_utils/fypp /bin/fypp

DBCSR_ver="2.6.0"

if [ -f dbcsr-${DBCSR_ver}.tar.gz ]; then
  echo "dbcsr-${DBCSR_ver}.tar.gz is found"
else
  wget -q "https://github.com/cp2k/dbcsr/archive/refs/tags/v${DBCSR_ver}.tar.gz" -O "dbcsr-${DBCSR_ver}.tar.gz"
fi

[ -d dbcsr-${DBCSR_ver} ] && rm -rf dbcsr-${DBCSR_ver}
tar xzf dbcsr-${DBCSR_ver}.tar.gz > dbcsr-build.log 2>&1
cd dbcsr-${DBCSR_ver}
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/opt/cp2k -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_SMM=blas .. > dbcsr-build.log 2>&1
make > dbcsr-build.log 2>&1 && make install > dbcsr-build.log 2>&1
cd ../..
mkdir build
cd build

if ! cmake \
  -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
  -Werror=dev \
  -DCP2K_USE_VORI=ON \
  -DCP2K_USE_COSMA=NO \
  -DCP2K_BLAS_VENDOR=OpenBLAS \
  -DCP2K_USE_SPGLIB=ON \
  -DCP2K_USE_LIBINT2=ON \
  -DCP2K_USE_LIBXC=ON \
  -DCP2K_USE_LIBTORCH=OFF \
  -DCP2K_USE_MPI=ON \
  -DCP2K_USE_MPI_F08=ON \
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
