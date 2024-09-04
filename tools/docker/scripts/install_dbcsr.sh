#!/bin/bash -e

DBCSR_ver="2.7.0"
DBCSR_sha256="25c367b49fb108c5230bcfb127f05fc16deff2bb467f437023dfa6045aff66f6"

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
