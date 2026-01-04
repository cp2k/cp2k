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

if [[ "${PROFILE}" =~ ^spack ]]; then
  eval "$(spack env activate myenv --sh)"
  # PyTorch's TorchConfig.cmake is buried in the Python site-packages directory
  Torch_DIR="$(dirname "$(find /opt/spack ! -type l -name TorchConfig.cmake | tail -n 1)")"
  export Torch_DIR
elif [[ "${PROFILE}" =~ ^toolchain ]]; then
  # shellcheck disable=SC1091
  source "${TOOLCHAIN_DIR}/install/setup"
elif [[ "${PROFILE}" =~ ^ubuntu ]] || [[ "${PROFILE}" =~ ^minimal ]]; then
  # Toolchain is only used to install DBCSR (and sometimes CMake).
  # shellcheck disable=SC1091
  source "${TOOLCHAIN_DIR}/install/setup"
fi

# Run CMake
[[ -d build ]] && rm -rf build
mkdir build
cd build || return 1

# TODO: Reconcile PROFILE/VERSION with CP2K_BUILD_OPTIONS in CMakeLists.txt
#
if [[ "${PROFILE}" == "spack" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_TBLITE=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "spack" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -DCP2K_USE_TBLITE=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "pdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="Debug" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "sdbg" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="Debug" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_MPI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "psmp" ]]; then
  # TODO Re-enable SIRIUS once performance regression of COSMA is fixed:
  # https://github.com/cp2k/cp2k/issues/4663
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_SIRIUS=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_MPI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_intel" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=OFF \
    -DCP2K_USE_MPI=ON \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_intel" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=OFF \
    -DCP2K_USE_MPI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_arm64" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_ACE=OFF \
    -DCP2K_USE_DEEPMD=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_cuda_"* ]] && [[ "${VERSION}" == "psmp" ]]; then
  # TODO: CP2K_USE_SPLA_GEMM_OFFLOADING should get enabled automatically.
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DDBCSR_DIR="${DBCSR_CUDA_ROOT}/lib/cmake/dbcsr" \
    -DSPLA_ROOT="${SPLA_ROOT}-cuda" \
    -DCP2K_WITH_GPU="${PROFILE:15}" \
    -DCP2K_USE_ACCEL=CUDA \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_SPLA_GEMM_OFFLOADING=ON \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_ACE=OFF \
    -DCP2K_USE_DEEPMD=OFF \
    -DCP2K_USE_PLUMED=OFF \
    -DCP2K_USE_LIBSMEAGOL=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_DFTD4=OFF \
    -DCP2K_USE_TBLITE=OFF \
    -DCP2K_USE_GREENX=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_VORI=OFF \
    -DCP2K_USE_TREXIO=OFF \
    -DCP2K_USE_MIMIC=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_hip_"* ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DDBCSR_DIR="${DBCSR_HIP_ROOT}/lib/cmake/dbcsr" \
    -DCP2K_WITH_GPU="${PROFILE:14}" \
    -DCP2K_USE_ACCEL=HIP \
    -DCP2K_USE_MPI=ON \
    -DCP2K_USE_LIBXC=ON \
    -DCP2K_USE_LIBINT2=ON \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_generic" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="Generic" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_DEEPMD=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_conventions" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="Conventions" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCMAKE_Fortran_COMPILER_LAUNCHER="redirect_gfortran_output.py" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_coverage" ]] && [[ "${VERSION}" == "psmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="Coverage" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "toolchain_asan" ]] && [[ "${VERSION}" == "psmp" ]]; then
  # TODO Re-enable GREENX. It currently leads to a heap-buffer-overflow
  # in `greenx_refine_pade()` at greenx_interface.F:80.
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="ASAN" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_DLAF=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_GREENX=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "ubuntu" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  # TODO Enable libxsmm and libspglib once we dropped GCC 8.
  # TODO fix spglib https://github.com/cp2k/cp2k/issues/3414
  # NOTE: libxc 5.2.3 is provided, CP2K requires libxc 7
  cmake \
    -GNinja \
    -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_USE_EVERYTHING=ON \
    -DCP2K_USE_ACE=OFF \
    -DCP2K_USE_DEEPMD=OFF \
    -DCP2K_USE_DFTD4=OFF \
    -DCP2K_USE_TBLITE=OFF \
    -DCP2K_USE_GREENX=OFF \
    -DCP2K_USE_LIBTORCH=OFF \
    -DCP2K_USE_LIBXC=OFF \
    -DCP2K_USE_MPI=OFF \
    -DCP2K_USE_PEXSI=OFF \
    -DCP2K_USE_SPGLIB=OFF \
    -DCP2K_USE_VORI=OFF \
    -DCP2K_USE_TREXIO=OFF \
    -DCP2K_USE_LIBXSMM=OFF \
    -Werror=dev \
    .. |& tee ./cmake.log
  CMAKE_EXIT_CODE=$?

elif [[ "${PROFILE}" == "minimal" ]] && [[ "${VERSION}" == "ssmp" ]]; then
  cmake \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
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

# EOF
