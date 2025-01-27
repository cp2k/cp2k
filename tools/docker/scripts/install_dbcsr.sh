#!/bin/bash -e

if (($# != 2)); then
  echo "Usage: install_dbcsr.sh <PROFILE> <VERSION>"
  exit 1
fi

PROFILE=$1
VERSION=$2

DBCSR_ver="2.8.0"
DBCSR_sha256="d55e4f052f28d1ed0faeaa07557241439243287a184d1fd27f875c8b9ca6bd96"

echo "==================== Installing DBCSR ===================="

if [[ "${PROFILE}" == "toolchain-cwd" ]]; then
  CWD=${PWD}
  TOOLCHAIN_DIR="${CWD}/tools/toolchain"
  INSTALL_PREFIX="${TOOLCHAIN_DIR}/install/dbcsr-${DBCSR_ver}/${VERSION}"
  DBCSR_DIR="${INSTALL_PREFIX}/lib/cmake/dbcsr"
  echo "export DBCSR_DIR=${INSTALL_PREFIX}/lib/cmake/dbcsr" >> "${TOOLCHAIN_DIR}/install/setup"
  if [[ -d "${DBCSR_DIR}" ]]; then
    echo "dbcsr-${DBCSR_ver} is already installed, skipping it."
    exit
  else
    [[ ! -d "${TOOLCHAIN_DIR}/build" ]] && mkdir -p "${TOOLCHAIN_DIR}/build"
    cd "${TOOLCHAIN_DIR}/build"
  fi
  PROFILE="toolchain"
else
  TOOLCHAIN_DIR="/opt/cp2k-toolchain"
  INSTALL_PREFIX="/opt/dbcsr"
fi

wget -q "https://github.com/cp2k/dbcsr/releases/download/v${DBCSR_ver}/dbcsr-${DBCSR_ver}.tar.gz"
echo "${DBCSR_sha256}  dbcsr-${DBCSR_ver}.tar.gz" | sha256sum --check > /dev/null

tar xzf dbcsr-${DBCSR_ver}.tar.gz
cd dbcsr-${DBCSR_ver}

mkdir build
cd build

if [[ "${PROFILE}" == "toolchain" ]]; then
  if [[ -f "${TOOLCHAIN_DIR}/install/setup" ]]; then
    # shellcheck disable=SC1091
    source "${TOOLCHAIN_DIR}/install/setup"
  else
    echo "ERROR: Toolchain setup file not found"
    exit 1
  fi
fi

if [[ "${VERSION}" == "ssmp" ]]; then
  USE_MPI="OFF"
elif [[ "${VERSION}" == "psmp" ]]; then
  USE_MPI="ON"
else
  echo "Unknown version: ${VERSION}."
  exit 1
fi

if ! cmake -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" -DCMAKE_INSTALL_LIBDIR=lib -DUSE_MPI=${USE_MPI} -DUSE_OPENMP=ON -DUSE_SMM=blas .. &> cmake.log; then
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

cd ../..
rm -rf "dbcsr-${DBCSR_ver}" "dbcsr-${DBCSR_ver}.tar.gz"

#EOF
