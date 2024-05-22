#!/bin/bash

# author: Ole Schuett

if (($# != 2)); then
  echo "Usage: build_cp2k_cmake.sh <PROFILE> <VERSION>"
  exit 1
fi

PROFILE=$1
VERSION=$2

echo "==================== Building CP2K ===================="

# Using Ninja because of https://gitlab.kitware.com/cmake/cmake/issues/18188

# Run CMake.
mkdir build
cd build || exit 1

if [[ "${PROFILE}" == "spack" ]]; then
  eval "$(spack env activate myenv --sh)"
elif [[ "${PROFILE}" == "toolchain" ]]; then
  # shellcheck disable=SC1091
  source /opt/cp2k-toolchain/install/setup
fi

# TODO: Reconcile PROFILE/VERSION with CP2K_BUILD_OPTIONS in CMakeLists.txt.
if [[ "${PROFILE}" == "spack" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_C_FLAGS="-fno-lto" \
    -DCMAKE_Fortran_FLAGS="-fno-lto" \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_USE_VORI=OFF \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_LIBINT2=OFF \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_MPI_F08=ON \
    -DCP2K_ENABLE_REGTESTS=ON \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_ENABLE_REGTESTS=ON \
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

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "sdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_DEBUG_MODE=ON \
    -DCP2K_ENABLE_REGTESTS=ON \
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
  # TODO Fix ELPA, COSMA, SIRIUS, QUIP, PEXSI, and Torch.
  # https://github.com/cp2k/cp2k/issues/3416
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_ENABLE_REGTESTS=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_VORI=ON \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_MPI_F08=ON \
    -DCP2K_USE_LIBXSMM=ON \
    -DCP2K_USE_SUPERLU=ON \
    -DCP2K_USE_PLUMED=ON \
    -DCP2K_USE_PEXSI=ON \
    -DCP2K_USE_SPLA=ON \
    -DCP2K_USE_METIS=ON \
    -DCP2K_USE_ELPA=OFF \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_SIRIUS=OFF \
    -DCP2K_USE_QUIP=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_DLAF=OFF \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "pdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_DEBUG_MODE=ON \
    -DCP2K_ENABLE_REGTESTS=ON \
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
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX=/opt/cp2k \
    -Werror=dev \
    -DCP2K_ENABLE_REGTESTS=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBXSMM=ON \
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
    -DCP2K_ENABLE_REGTESTS=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_FFTW3=ON \
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
    -DCP2K_ENABLE_REGTESTS=ON \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_LIBINT2=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_FFTW3=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
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
echo -en '\nCompiling cp2k...'
if ninja --verbose &> ninja.log; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 ninja.log
  mkdir -p /workspace/artifacts/
  cp ninja.log /workspace/artifacts/
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 1
fi

#EOF
