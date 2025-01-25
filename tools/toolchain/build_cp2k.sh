#!/bin/bash
#
# Build script for CP2K (GNU x86_64)
#
# Tested with: GNU 14.2.0, MPICH 4.2.3 and OpenMPI 5.0.6 on Linux clusters
#
# Last update: 24.01.2025
#
# shellcheck disable=SC1091

if [[ -n $1 ]]; then
  version="$1"
else
  version="ssmp"
fi

cwd=$PWD

[[ -z "${log_lines}" ]] && log_lines=1000
[[ -z "${maxtasks}" ]] && maxtasks=12
[[ -z "${mpi_implementation}" ]] && mpi_implementation="mpich"
[[ -z "${target_cpu}" ]] && target_cpu="native"

# Build CP2K toolchain

cd tools/toolchain || exit 1
rm -rf build

echo "Found GCC $(gcc -dumpfullversion) compiler"

if [[ "${version}" == "ssmp" ]]; then
  echo "Building for target CPU ${target_cpu} without MPI"
  ./install_cp2k_toolchain.sh -j "${maxtasks}" --mpi-mode=no --no-arch-files --log-lines="${log_lines}" --target-cpu="${target_cpu}" \
    --with-gcc --with-dftd4 --with-hdf5 --with-libtorch --with-trexio
  USE_MPI="OFF"
elif [[ "${version}" == "psmp" ]]; then
  echo "Building for target CPU ${target_cpu} using ${mpi_implementation}"
  ./install_cp2k_toolchain.sh --install-all -j "${maxtasks}" --no-arch-files --log-lines="${log_lines}" --target-cpu="${target_cpu}" \
    --with-gcc --with-"${mpi_implementation}" --with-libsmeagol
  USE_MPI="ON"
else
  echo "ERROR: Unknown version ${version} specified"
  exit 1
fi

if [[ -f install/setup ]]; then
  source install/setup
else
  echo "ERROR: Toolchain setup file not found"
  exit 1
fi

# DBCSR installation

echo "==================== Installing DBCSR ===================="

DBCSR_ver="2.8.0"
DBCSR_sha256="d55e4f052f28d1ed0faeaa07557241439243287a184d1fd27f875c8b9ca6bd96"

export DBCSR_DIR="${PWD}/install/dbcsr-${DBCSR_ver}/${version}/lib/cmake/dbcsr"

if [[ -d ${DBCSR_DIR} ]]; then

  echo "dbcsr-${DBCSR_ver} is already installed, skipping it."

else

  INSTALL_PREFIX="${PWD}/install/dbcsr-${DBCSR_ver}/${version}"

  cd build || exit 1

  wget -q "https://github.com/cp2k/dbcsr/releases/download/v${DBCSR_ver}/dbcsr-${DBCSR_ver}.tar.gz"
  echo "${DBCSR_sha256}  dbcsr-${DBCSR_ver}.tar.gz" | sha256sum --check > /dev/null

  tar -xzf dbcsr-${DBCSR_ver}.tar.gz
  cd dbcsr-${DBCSR_ver} || exit 1

  [[ -d build ]] && rm -rf build
  mkdir build
  cd build || exit 1

  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DUSE_MPI=${USE_MPI} \
    -DUSE_OPENMP=ON \
    -DUSE_SMM=blas \
    -Werror=dev \
    .. > cmake.log 2>&1 || (
    tail -n "${log_lines}" cmake.log
    exit 1
  )

  ninja --verbose > ninja.log 2>&1 || (
    tail -n "${log_lines}" ninja.log
    exit 1
  )

  ninja --verbose install > install.log 2>&1 || (
    tail -n "${log_lines}" install.log
    exit 1
  )

  cd ../.. || exit 1
  rm -rf "dbcsr-${DBCSR_ver}" "dbcsr-${DBCSR_ver}.tar.gz"
  cd .. || exit 1

fi

cd "${cwd}" || exit 1

echo "==================== Installing CP2K ===================="

source tools/toolchain/install/setup

INSTALL_PREFIX=${cwd}/install/${version}

if [[ -d "${INSTALL_PREFIX}" ]]; then

  echo "cp2k.${version} is already installed, skipping it."

else

  [[ -d build ]] && rm -rf build
  mkdir build
  cd build || exit 1

  if [[ "${version}" == "ssmp" ]]; then

    cmake \
      -GNinja \
      -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCP2K_BLAS_VENDOR=OpenBLAS \
      -DCP2K_USE_COSMA=OFF \
      -DCP2K_USE_DFTD4=ON \
      -DCP2K_USE_DLAF=OFF \
      -DCP2K_USE_FFTW3=ON \
      -DCP2K_USE_LIBINT2=ON \
      -DCP2K_USE_LIBTORCH=ON \
      -DCP2K_USE_LIBXC=ON \
      -DCP2K_USE_LIBXSMM=ON \
      -DCP2K_USE_MPI=OFF \
      -DCP2K_USE_MPI_F08=OFF \
      -DCP2K_USE_SPGLIB=ON \
      -DCP2K_USE_VORI=ON \
      -Werror=dev \
      .. > cmake.log 2>&1 || (
      tail -n "${log_lines}" cmake.log
      exit 1
    )

  elif [[ "${version}" == "psmp" ]]; then

    cmake \
      -GNinja \
      -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
      -DCMAKE_INSTALL_LIBDIR=lib \
      -DCP2K_BLAS_VENDOR=OpenBLAS \
      -DCP2K_USE_COSMA=ON \
      -DCP2K_USE_DFTD4=ON \
      -DCP2K_USE_DLAF=OFF \
      -DCP2K_USE_ELPA=ON \
      -DCP2K_USE_FFTW3=ON \
      -DCP2K_USE_LIBINT2=ON \
      -DCP2K_USE_LIBTORCH=ON \
      -DCP2K_USE_LIBXC=ON \
      -DCP2K_USE_LIBXSMM=ON \
      -DCP2K_USE_MPI=ON \
      -DCP2K_USE_MPI_F08=ON \
      -DCP2K_USE_PLUMED=ON \
      -DCP2K_USE_SIRIUS=OFF \
      -DCP2K_USE_SPGLIB=ON \
      -DCP2K_USE_SPLA=ON \
      -DCP2K_USE_VORI=ON \
      -Wno-error=dev \
      .. > cmake.log 2>&1 || (
      tail -n "${log_lines}" cmake.log
      exit 1
    )

  else

    echo "ERROR: Unknown version ${version} specified"
    exit 1

  fi

  ninja --verbose > ninja.log 2>&1 || (
    tail -n "${log_lines}" ninja.log
    exit 1
  )

  ninja --verbose install > install.log 2>&1 || (
    tail -n "${log_lines}" install.log
    exit 1
  )

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${INSTALL_PREFIX}/lib

fi

"${INSTALL_PREFIX}/bin/cp2k.${version}" -h -v

echo ""
echo " Run the following commands after a successful build"
echo ""
echo "   source ${cwd}/tools/toolchain/install/setup"
echo "   export LD_LIBRARY_PATH=${INSTALL_PREFIX}/lib:\$LD_LIBRARY_PATH"
echo "   export PATH=${INSTALL_PREFIX}/bin:\$PATH"
echo ""
echo " A CP2K regression test can be run with"
echo "   tests/do_regtest.py install/${version}/bin ${version}"

#EOF
