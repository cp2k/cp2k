#!/bin/bash

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: build_cp2k_cmake.sh <PROFILE>"
  exit 1
fi

PROFILE=$1

echo "==================== Building CP2K ===================="

# Using Ninja because of https://gitlab.kitware.com/cmake/cmake/issues/18188

# Run CMake.
mkdir build
cd build || exit 1

# TODO: Reconcile PROFILE with CP2K_BUILD_OPTIONS in CMakeLists.txt.
case ${PROFILE} in
  spack)
    eval "$(spack env activate myenv --sh)"
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
    ;;
  ubuntu)
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
    ;;
  minimal)
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
    ;;
  *)
    echo "Unknown profile ${PROFILE}."
    exit 1
    ;;
esac

if [ $CMAKE_EXIT_CODE -ne 0 ]; then
  echo -e "\nSummary: CMake failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

# Check for CMake warnings.
if grep -A5 'CMake Warning' ./cmake.log; then
  echo -e "\nSummary: Found CMake warnings."
  echo -e "Status: FAILED\n"
  exit 0
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
  exit 0
fi

#EOF
