#!/bin/bash -e

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

DBCSR_ver="2.6.0"
DBCSR_sha256="c67b02ff9abc7c1f529af446a9f01f3ef9e5b0574f220259128da8d5ca7e9dc6"

echo "==================== Installing DBCSR ===================="

wget -q "https://github.com/cp2k/dbcsr/releases/download/v${DBCSR_ver}/dbcsr-${DBCSR_ver}.tar.gz"
echo "${DBCSR_sha256}  dbcsr-${DBCSR_ver}.tar.gz" | sha256sum --check > /dev/null

tar xzf dbcsr-${DBCSR_ver}.tar.gz
cd dbcsr-${DBCSR_ver}

mkdir build
cd build

if ! cmake -DCMAKE_INSTALL_PREFIX=/opt/dbcsr -DUSE_MPI=ON -DUSE_OPENMP=ON -DUSE_SMM=blas .. &> cmake.log; then
  cat cmake.log
  exit 1
fi

if ! make -j &> make.log; then
  cat make.log
  exit 1
fi

if ! make install VERBOSE=1 &> install.log; then
  cat install.log
  exit 1
fi

cd ..
rm -rf build "dbcsr-${DBCSR_ver}" "dbcsr-${DBCSR_ver}.tar.gz"

#EOF
