# Common variables used by the installation scripts

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all
# shellcheck shell=bash

# directories and files used by the installer
ROOTDIR=${ROOTDIR:-"$(pwd -P)"}
SCRIPTDIR=${SCRIPTDIR:-"${ROOTDIR}/scripts"}
INSTALLDIR=${INSTALLDIR:-"${ROOTDIR}/install"}
BUILDDIR=${BUILDDIR:-"${ROOTDIR}/build"}
SETUPFILE=${SETUPFILE:-"${INSTALLDIR}/setup"}
VERSION_FILE=${VERSION_FILE:-"${SCRIPTDIR}/VERSION"}

# Normalize the system architecture reported by uname across supported platforms.
SYSTEM_ARCH=${SYSTEM_ARCH:-"$(uname -m)"}
case "${SYSTEM_ARCH}" in
  x86_64 | amd64 | x86-64)
    SYSTEM_ARCH="x86_64"
    ;;
  aarch64 | arm64)
    SYSTEM_ARCH="arm64"
    ;;
  *)
    # Keep other architectures as reported by uname.
    ;;
esac
export SYSTEM_ARCH

# search paths
SYS_INCLUDE_PATH=${SYS_INCLUDE_PATH:-'/usr/local/include:/usr/include'}
SYS_LIB_PATH=${SYS_LIB_PATH:-"/usr/local/lib64:/usr/local/lib:/usr/lib64:/usr/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/aarch64-linux-gnu:/lib64:/lib"}
INCLUDE_PATHS=${INCLUDE_PATHS:-"CPATH SYS_INCLUDE_PATH"}
LIB_PATHS=${LIB_PATHS:-'LD_LIBRARY_PATH LIBRARY_PATH LD_RUN_PATH SYS_LIB_PATH'}

# mode flags
ENABLE_OMP=${ENABLE_OMP:-"__TRUE__"}
ENABLE_CUDA=${ENABLE_CUDA:-"__FALSE__"}
ENABLE_GAUXC_CUTLASS=${ENABLE_GAUXC_CUTLASS:-"__FALSE__"}
ENABLE_HIP=${ENABLE_HIP:-"__FALSE__"}
ENABLE_OPENCL=${ENABLE_OPENCL:-"__FALSE__"}
ENABLE_CRAY=${ENABLE_CRAY:-"__FALSE__"}
MPI_MODE=${MPI_MODE:-openmpi}
MATH_MODE=${MATH_MODE:-openblas}

export NVCC=${NVCC:-nvcc}
