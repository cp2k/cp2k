#!/bin/bash

# Authors: Ole Schuett, Matthias Krack

if [[ "${0}" == "${BASH_SOURCE[0]}" ]]; then
  echo "ERROR: Script ${0##*/} must be sourced"
  echo "Usage: ${BASH_SOURCE##*/} <PROFILE> <VERSION>"
  exit 1
fi

if (($# != 2)); then
  echo "ERROR: Script ${BASH_SOURCE##*/} expects exactly two arguments"
  echo "Usage: ${BASH_SOURCE##*/} <PROFILE> <VERSION>"
  return 1
fi

PROFILE=$1
VERSION=$2

[[ -z "${TOOLCHAIN_DIR}" ]] && TOOLCHAIN_DIR="/opt/cp2k-toolchain"
[[ -z "${INSTALL_PREFIX}" ]] && INSTALL_PREFIX="/opt/cp2k"

# Using Ninja because of https://gitlab.kitware.com/cmake/cmake/issues/18188

if [[ "${PROFILE}" == "spack" ]]; then
  eval "$(spack env activate myenv --sh)"
elif [[ "${PROFILE}" == "toolchain" ]]; then
  # shellcheck disable=SC1091
  source "${TOOLCHAIN_DIR}/install/setup"
fi

# Run CMake
[[ -d build ]] && rm -rf build
mkdir build
cd build || return 1

# TODO: Reconcile PROFILE/VERSION with CP2K_BUILD_OPTIONS in CMakeLists.txt
if [[ "${PROFILE}" == "spack" ]] && [[ "${VERSION}" == "psmp" ]]; then
  # PyTorch's TorchConfig.cmake is buried in the Python site-packages directory
  export Torch_DIR="/opt/spack/var/spack/environments/myenv/spack-env/view/lib/python3.11/site-packages/torch/share/cmake/Torch"
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_C_FLAGS="-fno-lto" \
    -DCMAKE_Fortran_FLAGS="-fno-lto" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR="OpenBLAS" \
    -DCP2K_SCALAPACK_VENDOR="GENERIC" \
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
    -DCP2K_USE_LIBVDWXC=ON \
    -DCP2K_USE_GRPP=OFF \
    -DCP2K_USE_TREXIO=ON \
    -DCP2K_USE_LIBTORCH=ON \
    -DCP2K_USE_DLAF=ON \
    -DCP2K_USE_DFTD4=ON \
    -DCP2K_USE_LIBSMEAGOL=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DEEPMD=ON \
    -DCP2K_USE_DFTD4=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_GRPP=OFF \
    -DCP2K_USE_HDF5=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBTORCH=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBXSMM=ON \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_TREXIO=ON \
    -DCP2K_USE_VORI=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "sdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_DEBUG_MODE=ON \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_VORI=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_LIBDIR=lib \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=ON \
    -DCP2K_USE_DEEPMD=ON \
    -DCP2K_USE_DFTD4=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_ELPA=ON \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_GRPP=OFF \
    -DCP2K_USE_HDF5=ON \
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
    -DCP2K_USE_TREXIO=ON \
    -DCP2K_USE_VORI=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "pdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_DEBUG_MODE=ON \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_MPI_F08=ON \
    -DCP2K_USE_SPGLIB=ON \
    -DCP2K_USE_VORI=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "ubuntu" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  # TODO fix spglib https://github.com/cp2k/cp2k/issues/3414
  # NOTE: libxc 5.2.3 is provided, CP2K requires libxc 7
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_LIBXSMM=ON \
    -DCP2K_USE_HDF5=ON \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_VORI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "ubuntu_i386" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  # TODO fix spglib https://github.com/cp2k/cp2k/issues/3414
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_FFTW3=ON \
    -DCP2K_USE_LIBINT2=ON \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_VORI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
elif [[ "${PROFILE}" == "minimal" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_BLAS_VENDOR=OpenBLAS \
    -DCP2K_USE_COSMA=OFF \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_FFTW3=OFF \
    -DCP2K_USE_LIBINT2=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_MPI_F08=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_VORI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?
else
  echo "Unknown combination of PROFILE=\"${PROFILE}\" and VERSION=\"${VERSION}\""
  return 1
fi

if ((CMAKE_EXIT_CODE != 0)); then
  echo -e "\nSummary: CMake failed"
  echo -e "Status: FAILED\n"
  return 1
fi

# Check for CMake warnings
if grep -A5 'CMake Warning' ./cmake.log; then
  echo -e "\nSummary: Found CMake warnings"
  echo -e "Status: FAILED\n"
  return 1
fi

#EOF
