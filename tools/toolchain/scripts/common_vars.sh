# Common variables used by the installation scripts

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237
# shellcheck shell=bash

# directories and files used by the installer
ROOTDIR=${ROOTDIR:-"$(pwd -P)"}
SCRIPTDIR=${SCRIPTDIR:-"${ROOTDIR}/scripts"}
INSTALLDIR=${INSTALLDIR:-"${ROOTDIR}/install"}
BUILDDIR=${BUILDDIR:-"${ROOTDIR}/build"}
SETUPFILE=${SETUPFILE:-"${INSTALLDIR}/setup"}
ARCH_FILE_TEMPLATE=${ARCH_FILE_TEMPLATE:-"${SCRIPTDIR}/arch_base.tmpl"}
VERSION_FILE=${VERSION_FILE:-"${SCRIPTDIR}/VERSION"}

# downloader flags, used for downloading tarballs, see download_pkg macro in tool_kit.sh
DOWNLOADER_FLAGS="${DOWNLOADER_FLAGS:-}"

# system arch gotten from OpenBLAS prebuild
OPENBLAS_ARCH=${OPENBLAS_ARCH:-"x86_64"}
OPENBLAS_LIBCORE=${OPENBLAS_LIBCORE:-''}

# search paths
SYS_INCLUDE_PATH=${SYS_INCLUDE_PATH:-'/usr/local/include:/usr/include'}
SYS_LIB_PATH=${SYS_LIB_PATHS:-'/usr/local/lib64:/usr/local/lib:/usr/lib64:/usr/lib:/lib64:/lib'}
INCLUDE_PATHS=${INCLUDE_PATHS:-"CPATH SYS_INCLUDE_PATH"}
LIB_PATHS=${LIB_PATHS:-'LD_LIBRARY_PATH LIBRARY_PATH LD_RUN_PATH SYS_LIB_PATH'}

# mode flags
ENABLE_OMP=${ENABLE_OMP:-"__TRUE__"}
ENABLE_CUDA=${ENABLE_CUDA:-"__FALSE__"}
ENABLE_CRAY=${ENABLE_CRAY:-"__FALSE__"}
MPI_MODE=${MPI_MODE:-openmpi}
MATH_MODE=${MATH_MODE:-openblas}

# compiler flags
export CC=${CC:-gcc}
export FC=${FC:-gfortran}
export F77=${F77:-gfortran}
export F90=${F90:-gfortran}
export CXX=${CXX:-g++}
export MPICC=${MPICC:-mpicc}
export MPIFC=${MPIFC:-mpif90}
export MPIF77=${MPIF77:-mpif77}
export MPIF90=${MPIF90:-mpif90}
export MPICXX=${MPICXX:-mpic++}
export NVCC=${NVCC:-nvcc}
