#!/bin/bash

# authors: Ole Schuett, Matthias Krack

if (($# != 2)); then
  echo "Usage: build_cp2k_cmake.sh <PROFILE> <VERSION>"
  exit 1
fi

PROFILE=$1
VERSION=$2

[[ -z "${log_lines}" ]] && log_lines="100"
[[ -z "${maxtasks}" ]] && maxtasks="12"
[[ -z "${mpi_implementation}" ]] && mpi_implementation="mpich"
[[ -z "${target_cpu}" ]] && target_cpu="native"

if [[ "${PROFILE}" == "toolchain-cwd" ]]; then
  CWD="${PWD}"
  TOOLCHAIN_DIR="${CWD}/tools/toolchain"
  INSTALL_PREFIX="${CWD}/install/${VERSION}"
  # Build or update CP2K toolchain if needed
  echo "==================== Building CP2K toolchain =========="
  cd "${TOOLCHAIN_DIR}" || exit 1
  rm -rf build
  echo "Found GCC $(gcc -dumpfullversion) compiler"
  if [[ "${VERSION}" == "ssmp" ]]; then
    echo "Building for target CPU ${target_cpu} without MPI"
    ./install_cp2k_toolchain.sh -j "${maxtasks}" --mpi-mode=no --no-arch-files --log-lines="${log_lines}" --target-cpu="${target_cpu}" \
      --with-gcc --with-dftd4 --with-hdf5 --with-libtorch --with-trexio
  elif [[ "${VERSION}" == "psmp" ]]; then
    echo "Building for target CPU ${target_cpu} using ${mpi_implementation}"
    ./install_cp2k_toolchain.sh --install-all -j "${maxtasks}" --no-arch-files --log-lines="${log_lines}" --target-cpu="${target_cpu}" \
      --with-gcc --with-"${mpi_implementation}" --with-libsmeagol
  else
    echo "ERROR: Unknown version ${VERSION} specified"
    exit 1
  fi
  cd "${CWD}" || exit 1
  # Install DBCSR if needed
  "${CWD}/tools/docker/scripts/install_dbcsr.sh" "${PROFILE}" "${VERSION}"
  PROFILE="toolchain"
else
  TOOLCHAIN_DIR="/opt/cp2k-toolchain"
  INSTALL_PREFIX="/opt/cp2k"
fi

echo "==================== Building CP2K ===================="

# Using Ninja because of https://gitlab.kitware.com/cmake/cmake/issues/18188

# Run CMake.
[[ -d build ]] && rm -rf build
mkdir build
cd build || exit 1

if [[ "${PROFILE}" == "spack" ]]; then
  eval "$(spack env activate myenv --sh)"
elif [[ "${PROFILE}" == "toolchain" ]]; then
  # shellcheck disable=SC1091
  source "${TOOLCHAIN_DIR}/install/setup"
fi

# TODO: Reconcile PROFILE/VERSION with CP2K_BUILD_OPTIONS in CMakeLists.txt.
if [[ "${PROFILE}" == "spack" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_C_FLAGS="-fno-lto" \
    -DCMAKE_Fortran_FLAGS="-fno-lto" \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_HDF5=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_VORI=ON \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_MPI_F08=ON \
    -DCP2K_USE_LIBXSMM=ON \
    -DCP2K_USE_PLUMED=ON \
    -DCP2K_USE_SPLA=ON \
    -DCP2K_USE_ELPA=ON \
    -DCP2K_USE_COSMA=ON \
    -DCP2K_USE_SIRIUS=ON \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_DLAF=ON \
    -DCP2K_USE_DFTD4=ON \
    -DCP2K_USE_LIBSMEAGOL=ON \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DEEPMD=ON \
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
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "sdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_DEBUG_MODE=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_VORI=ON \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=ON \
    -DCP2K_USE_DEEPMD=ON \
    -DCP2K_USE_DFTD4=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_ELPA=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBSMEAGOL=ON \
    -DCP2K_USE_LIBTORCH=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBXSMM=ON \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_MPI_F08=ON \
    -DCP2K_USE_PLUMED=ON \
    -DCP2K_USE_SIRIUS=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_SPLA=ON \
    -DCP2K_USE_VORI=ON \
    -DCP2K_USE_HDF5=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "pdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_DEBUG_MODE=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_MPI_F08=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_VORI=ON \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "ubuntu" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  # TODO fix spglib https://github.com/cp2k/cp2k/issues/3414
  # NOTE: libxc 5.2.3 is provided, CP2K requires libxc 7
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBXSMM=ON \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_HDF5=ON \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_VORI=OFF \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "ubuntu_i386" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  # TODO fix spglib https://github.com/cp2k/cp2k/issues/3414
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_VORI=OFF \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "minimal" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_FFTW3=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_VORI=OFF \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

else
  echo "Unknown combination profile: ${PROFILE} version: ${VERSION}."
  exit 1
fi

if [ $CMAKE_EXIT_CODE -ne 0 ]; then
  echo -e "\nSummary: CMake failed."
  echo -e "Status: FAILED\n"
  exit 1
fi

# Check for CMake warnings.
if grep -A5 'CMake Warning' ./cmake.log; then
  echo -e "\nSummary: Found CMake warnings."
  echo -e "Status: FAILED\n"
  exit 1
fi

# Compile CP2K.
echo -en '\nCompiling CP2K ... '
if ninja --verbose &> ninja.log; then
  echo "done."
  # Install CP2K.
  echo -en '\nInstalling CP2K ... '
  if ninja --verbose install &> install.log; then
    echo "done."
    echo ""
    echo " Run the following commands after a successful CP2K build"
    echo ""
    echo "   source ${CWD}/tools/toolchain/install/setup"
    echo "   export LD_LIBRARY_PATH=${INSTALL_PREFIX}/lib:\$LD_LIBRARY_PATH"
    echo "   export PATH=${INSTALL_PREFIX}/bin:\$PATH"
    echo ""
    echo " A CP2K regression test can be run with"
    echo "   tests/do_regtest.py install/${VERSION}/bin ${VERSION}"
  else
    echo -e "failed.\n\n"
    tail -n "${log_lines}" install.log
    exit 1
  fi
else
  echo -e "failed.\n\n"
  tail -n "${log_lines}" ninja.log
  mkdir -p /workspace/artifacts/
  cp ninja.log /workspace/artifacts/
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 1
fi

#EOF
