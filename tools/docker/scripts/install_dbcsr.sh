#!/bin/bash -e

DBCSR_ver="2.8.0"
DBCSR_sha256="d55e4f052f28d1ed0faeaa07557241439243287a184d1fd27f875c8b9ca6bd96"

echo "==================== Installing DBCSR ===================="

wget -q "https://github.com/cp2k/dbcsr/releases/download/v${DBCSR_ver}/dbcsr-${DBCSR_ver}.tar.gz"
echo "${DBCSR_sha256}  dbcsr-${DBCSR_ver}.tar.gz" | sha256sum --check > /dev/null

tar xzf dbcsr-${DBCSR_ver}.tar.gz
cd dbcsr-${DBCSR_ver}

mkdir build
cd build

if ! cmake -DCMAKE_INSTALL_PREFIX=/opt/dbcsr -DUSE_MPI=OFF -DUSE_OPENMP=ON -DUSE_SMM=blas .. &> cmake.log; then
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
