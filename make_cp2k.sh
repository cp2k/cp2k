#!/usr/bin/env bash

# Purpose: Build CP2K using Spack and CMake locally within the folder CP2K_ROOT
#          which defaults to the current working directory. This should be a
#          "cp2k/" folder containing the CP2K source tree.
#
#          The script can either be sourced with
#
#          > source ./make_cp2k.sh
#
#          or run in a subshell with
#
#          > ./make_cp2k.sh
#
#          The latter, running in a subshell, is recommended.
#
#          The flags -h or --help print the available options.
#
#          The first run will take longer as it will build all CP2K dependencies
#          with Spack. The Spack installation is kept fully local in the subfolder
#          cp2k/spack which corresponds to the tools/toolchain/install folder
#          created by the CP2K toolchain.
#
#          Subsequent runs of the script will use the software stack from that
#          cp2k/spack/ folder.
#
#          A rebuild of all CP2K dependencies can be enforced simply by removing
#          or renaming the folder cp2k/spack. The latter allows for keeping different
#          software stacks (see also -bd and -bd_only flags).
#
#          It is recommended to install podman to take advantage of a local cache.
#          This will accelerate the (re)build of the CP2K dependencies with Spack
#          significantly.
#
#          After the CP2K dependencies are built with Spack, CP2K itself is built
#          and installed using CMake in the subfolders cp2k/build and cp2k/install,
#          respectively.
#
#          Subsequent runs of the script will use the CMake configuration in the
#          subfolder cp2k/build. A rebuild of CP2K from scratch can be enforced
#          by removing or renaming that subfolder.
#
#          A CP2K regression run can be launched automatically by adding the flag
#          -t "" (or --test ""). This flag expects a string with the TESTOPTS, e.g.
#
#          > ./make_cp2k.sh -t "--maxtasks 8 --restrictdir QS/regtest-gpw-1"
#
#          Alternatively, the script cp2k/install/run_tests can be launched after
#          a successful CP2K build.

# Authors: Matthias Krack (MK)

# Version: 1.5

# History: - Creation (19.12.2025, MK)
#          - Version 0.1: First working version (09.01.2026, MK)
#          - Version 0.2: Add more flags and checks (19.01.2026, MK)
#          - Version 0.3: Add no_externals flag and perform more checks (21.01.2026, MK)
#          - Version 0.4: Improve error handling and provide more hints (22.01.2026, MK)
#          - Version 0.5: Adapt script for use within a container (24.01.2026, MK)
#          - Version 0.6: Add MPI flag and revise flag parsing (27.01.2026, MK)
#          - Version 0.7: Fix container detection (28.01.2026, MK)
#          - Version 0.8: Add --build_deps_only flag (29.01.2026, MK)
#          - Version 0.9: Add --disable_local_cache flag (30.01.2026, MK)
#          - Version 1.0: Add Cray specific configuration (01.02.2026, MK)
#          - Version 1.1: Allow for selecting the GCC version (02.02.2026, MK)
#          - Version 1.2: Add option for static build (05.02.2026, MK)
#          - Version 1.3: Add CUDA GPU support (10.02.2026, MK)
#          - Version 1.4: Drop download of spack-packages (12.02.2026, MK)
#          - Version 1.5: Add flags to enable/disable features selectively (15.02.2026, MK)

# Facilitate the deugging of this script
set -uo pipefail

# Retrieve script name
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"

# Check if the script is sourced or run in a subshell
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
  echo "${SCRIPT_NAME}: Running script in sourcing mode"
  echo "${SCRIPT_NAME}: All changes from this script will become effective within the current shell"
  EXIT_CMD="return"
else
  echo "${SCRIPT_NAME}: Running script in subshell mode"
  echo "${SCRIPT_NAME}: Changes from this script in the shell environment are lost after completion"
  EXIT_CMD="exit"
fi

# Help finding ldconfig
if ! command -v ldconfig &> /dev/null; then
  PATH="${PATH}:/sbin"
fi

# Check if all mandatory packages are installed in the environment
for package in awk bzip2 find g++ gcc gfortran git gzip ldconfig make patch python3 tar wget xz; do
  if ! command -v "${package}" &> /dev/null; then
    echo "ERROR: The package \"${package}\" is mandatory to build CP2K with Spack/CMake"
    echo "       Install the missing package and re-run the script"
    ${EXIT_CMD} 1
  fi
done

# Check platform
if [[ "$(uname -s)" == "Darwin" ]]; then
  echo "ERROR: A build for Apple macOS (Darwin, Homebrew) is not supported yet,"
  echo "       but CP2K can be built within a container, e.g. using Ubuntu"
  ${EXIT_CMD} 1
fi

# Check bash version
if ((BASH_VERSINFO < 4)); then
  echo "ERROR: The employed bash version ${BASH_VERSION} is too old (from 2004)"
  echo "       Install a newer bash version"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    echo "HINT:  Use the bash command provided by Homebrew with Apple"
    echo "       /opt/homebrew/bin should precede /bin/bash in the PATH"
    echo "       In sourcing mode, run /opt/homebrew/bin/bash first before sourcing this script"
  fi
  ${EXIT_CMD} 1
elif ((BASH_VERSINFO < 5)); then
  echo "WARNING: The employed bash version ${BASH_VERSION} is quite old (from 2009)"
else
  echo "INFO: Using bash version ${BASH_VERSION}"
fi

# Check if the python3 version is new enough for spack
if ! python3 -c 'import sys; sys.exit(not(sys.version_info >= (3, 10)))'; then
  echo "ERROR: Python version is NOT >= 3.10 (needed for Spack)"
  echo "       Found only $(python3 -V)"
  ${EXIT_CMD} 1
else
  echo "INFO: Found $(python3 -V)"
fi

# Default values
BUILD_DEPS="if_needed"
BUILD_DEPS_ONLY="no"
BUILD_TYPE="${BUILD_TYPE:-Release}"
CMAKE_FEATURE_FLAG_ALL="-DCP2K_USE_EVERYTHING=ON" # all features are activated by default
CMAKE_FEATURE_FLAGS="-DCP2K_BLAS_VENDOR=OpenBLAS" # LAPACK/BLAS from OpenBLAS by default
CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_FFTW3=ON"       # FFTW3 is always activated unless explicitly disabled
CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_DLAF=OFF"       # DLAF is deactivated by default
CMAKE_FEATURE_FLAG_MPI="-DCP2K_USE_MPI=ON"        # MPI is switched on by default
CMAKE_FEATURE_FLAGS_GPU="-DCP2K_USE_SPLA_GEMM_OFFLOADING=ON"
CRAY="no"
CUDA_SM_CODE=0
DISABLE_LOCAL_CACHE="no"
GCC_VERSION="auto"
GPU_MODEL="none"
HAS_PODMAN="no"
HELP="no"
INSTALL_MESSAGE="NEVER"
MPI_MODE="mpich"
if command -v nproc &> /dev/null; then
  MAX_PROCS=$(nproc)
  NUM_PROCS=${NUM_PROCS:-${MAX_PROCS}}
else
  MAX_PROCS=-1
  NUM_PROCS=${NUM_PROCS:-8}
fi
NVCC_VERSION=0
REBUILD_CP2K="no"
RUN_TEST="no"
SED_PATTERN_LIST=""
TESTOPTS=""
USE_EXTERNALS="no"
VERBOSE=0
VERBOSE_FLAG="--quiet"
VERBOSE_MAKEFILE="OFF"

export CP2K_ENV="cp2k_env"
export CP2K_ROOT=${CP2K_ROOT:-${PWD}}
export CP2K_VERSION="${CP2K_VERSION:-psmp}"
export INSTALL_PREFIX="${INSTALL_PREFIX:-${CP2K_ROOT}/install}"

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -bd | --build_deps | --build_dependencies)
      BUILD_DEPS="always"
      shift 1
      ;;
    -bd_only | --build_deps_only | --build_dependencies_only)
      BUILD_DEPS="always"
      BUILD_DEPS_ONLY="yes"
      shift 1
      ;;
    -bt | --build_type)
      BUILD_TYPE="${2}"
      shift 2
      ;;
    -cray)
      CRAY="yes"
      shift 1
      ;;
    -cv | --cp2k_version)
      if (($# > 1)); then
        case "${2,,}" in
          psmp | ssmp | ssmp-static)
            CP2K_VERSION="${2,,}"
            ;;
          *)
            echo "ERROR: Invalid CP2K version \"${2}\" specified (choose psmp, ssmp, or ssmp-static)"
            ${EXIT_CMD} 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\" (choose psmp, ssmp, or ssmp-static)"
        ${EXIT_CMD} 1
      fi
      # Disable MPI for a serial CP2K binary
      if [[ "${CP2K_VERSION}" == "ssmp"* ]]; then
        MPI_MODE="no"
        CMAKE_FEATURE_FLAG_MPI="-DCP2K_USE_MPI=OFF"
      fi
      if [[ "${CP2K_VERSION}" == "ssmp-static" ]]; then
        CMAKE_FEATURE_FLAG_ALL="-DCP2K_USE_EVERYTHING=ON"
        for package in dftd4 libint2 libxc libxsmm spglib vori tblite; do
          CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_${package^^}=ON"
        done
        for package in ace deepmd greenx hdf5 libtorch pexsi trexio; do
          CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_${package^^}=OFF"
        done
      fi
      shift 2
      ;;
    -df | --disable | --disable_feature | -ef | --enable | --enable_feature)
      if (($# > 1)); then
        case "${1}" in
          -df | --disable | --disable_feature)
            ON_OFF="OFF"
            SUBST="s/^ /#/'"
            ;;
          -ef | --enable | --enable_feature)
            ON_OFF="ON"
            SUBST="s/#/ /'"
            ;;
        esac
        case "${2,,}" in
          all)
            CMAKE_FEATURE_FLAG_ALL="-DCP2K_USE_EVERYTHING=${ON_OFF}"
            for package in adios2 cosma deepmdkit dftd4 dla-future dla-future-fortran \
              elpa greenx hdf5 libfabric libint libvdwxc libsmeagol libvori libxc \
              libxsmm mimic-mcl openpmd-api pace pexsi plumed py-torch sirius spfft \
              spglib spla tblite trexio; do
              SED_PATTERN_LIST+=" -e '/\s*-\s+\"${package}@/ ${SUBST}"
            done
            if [[ "${ON_OFF}" == "OFF" ]]; then
              SED_PATTERN_LIST+=" -e '/\s*-\s+\"smm=libxsmm\"/ s/libxsmm/blas/'"
            fi
            ;;
          ace | cosma | deepmd | dftd4 | dlaf | elpa | fftw3 | greenx | hdf5 | libint2 | \
            libsmeagol | libtorch | libxc | libxsmm | mimic | openpmd | pexsi | plumed | \
            spglib | tblite | trexio | vori)
            CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_${2^^}=${ON_OFF}"
            # Translate package selection to sed pattern
            case "${2,,}" in
              ace)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"p${2,,}@/ ${SUBST}"
                ;;
              cosma | dftd4 | elpa | greenx | hdf5 | libsmeagol | libxc | pexsi | plumed | spglib | tblite | trexio)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}@/ ${SUBST}"
                ;;
              deepmd | libtorch)
                CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_DEEPMD=${ON_OFF}"
                CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_LIBTORCH=${ON_OFF}"
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}kit@/ ${SUBST}"
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"py-torch@/ ${SUBST}"
                ;;
              dlaf)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"dla-future.*@/ ${SUBST}"
                ;;
              fftw3)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"fftw@/ ${SUBST}"
                ;;
              libint2)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"libint@/ ${SUBST}"
                ;;
              libxsmm)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}@/ ${SUBST}"
                if [[ "${ON_OFF}" == "OFF" ]]; then
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"smm=${2,,}\"/ s/${2,,}/blas/'"
                fi
                ;;
              mimic)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"mimic-mcl@/ ${SUBST}"
                ;;
              openpmd | adios2)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"adios2@/ ${SUBST}"
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"openpmd-api@/ ${SUBST}"
                ;;
              vori)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"lib${2,,}@/ ${SUBST}"
                ;;
            esac
            ;;
          libvdwxc | spfft | spla | sirius)
            echo "WARNING: You have enabled or disabled one of the packages libvdwxc, spfft, spla, sirius"
            echo "         which means that the other packages will also be enabled or disabled"
            CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_LIBVDWXC=${ON_OFF} -DCP2K_USE_SpFFT=${ON_OFF}"
            CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_SPLA=${ON_OFF} -DCP2K_USE_SIRIUS=${ON_OFF}"
            SED_PATTERN_LIST+=" -e '/\s*-\s+\"libvdwxc@/ ${SUBST}"
            SED_PATTERN_LIST+=" -e '/\s*-\s+\"spfft@/ ${SUBST}"
            SED_PATTERN_LIST+=" -e '/\s*-\s+\"spla@/ ${SUBST}"
            SED_PATTERN_LIST+=" -e '/\s*-\s+\"sirius@/ ${SUBST}"
            ;;
          cray_pm_accel_energy | cusolver_mp | spla_gemm_offloading | unified_memory)
            CMAKE_FEATURE_FLAGS_GPU+=" -DCP2K_USE_${2^^}=${ON_OFF}"
            ;;
          dbm_gpu | elpa_gpu | grid_gpu | pw_gpu)
            CMAKE_FEATURE_FLAGS_GPU+=" -DCP2K_ENABLE_${2^^}=${ON_OFF}"
            ;;
          none)
            # do nothing
            ;;
          *)
            echo "ERROR: Unknown CP2K feature \"${2}\" specified"
            ${EXIT_CMD} 1
            ;;
        esac
      else
        echo "ERROR: No feature found for flag \"${1}\""
        ${EXIT_CMD} 1
      fi
      shift 2
      ;;
    -dlc | --disable_local_cache)
      DISABLE_LOCAL_CACHE="yes"
      shift 1
      ;;
    -gv | --gcc_version)
      if (($# > 1)); then
        case "${2}" in
          1[0-6])
            GCC_VERSION="${2}"
            ;;
          *)
            echo "ERROR: Invalid GCC version \"${2}\" specified (choose from 10 to 16)"
            ${EXIT_CMD} 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\" (choose a GCC version)"
        ${EXIT_CMD} 1
      fi
      shift 2
      ;;
    -gm | -gpu | --gpu_model)
      if (($# > 1)); then
        case "${2^^}" in
          P100 | V100 | T400 | A100 | A40 | H100 | H200 | GH200)
            GPU_MODEL="${2^^}"
            case "${GPU_MODEL}" in
              P100)
                CUDA_SM_CODE=60
                ;;
              V100)
                CUDA_SM_CODE=70
                ;;
              T400)
                CUDA_SM_CODE=75
                ;;
              A100)
                CUDA_SM_CODE=80
                ;;
              A40)
                CUDA_SM_CODE=86
                ;;
              H100 | H200 | GH200)
                CUDA_SM_CODE=90
                ;;
            esac
            ;;
          60 | 70 | 75 | 80 | 86 | 87 | 89 | 90 | 120 | 121)
            CUDA_SM_CODE=${2}
            ;;
          NONE)
            GPU_MODEL="${2,,}"
            ;;
          *)
            echo -e "\nERROR: Unknown GPU model \"${2}\" specified (choose <CUDA SM code>, P100, V100, T400, A100, A40, H100, H200, GH200 or none)\n"
            ${EXIT_CMD} 1
            ;;
        esac
        if ((CUDA_SM_CODE > 0)); then
          # Currently needed
          USE_EXTERNALS="yes"
          echo "INFO: The use of externals (-ue flag) is currently enforced with CUDA"
        fi
      else
        echo -e "\nERROR: No argument found for flag \"${1}\" (choose <CUDA SM code>, P100, V100, T400, A100, A40, H100, H200, GH200 or none)\n"
        ${EXIT_CMD} 1
      fi
      shift 2
      ;;
    -h | --help)
      HELP="yes"
      shift 1
      ;;
    -ip | --install_path | --install_prefix)
      if (($# > 1)); then
        INSTALL_PREFIX="${2}"
      else
        echo "ERROR: No install path argument found for flag \"${1}\""
        ${EXIT_CMD} 1
      fi
      shift 2
      ;;
    -j)
      if (($# > 1)); then
        case "${2}" in
          -*)
            shift 1
            ;;
          [0-9]*)
            NUM_PROCS="${2}"
            shift 2
            ;;
          *)
            echo "ERROR: The -j flag can only be followed by an integer number, found \"${2}\""
            ${EXIT_CMD} 1
            ;;
        esac
      else
        shift 1
      fi
      ;;
    -j[0-9]*)
      NUM_PROCS="${1#-j}"
      shift 1
      ;;
    -mpi | --mpi_mode)
      if (($# > 1)); then
        case "${2,,}" in
          mpich | openmpi)
            CMAKE_FEATURE_FLAG_MPI="-DCP2K_USE_MPI=ON"
            MPI_MODE="${2,,}"
            ;;
          no | none | off)
            CMAKE_FEATURE_FLAG_MPI="-DCP2K_USE_MPI=OFF"
            MPI_MODE="no"
            ;;
          *)
            echo "ERROR: Unknown MPI mode \"${2}\" specified (choose mpich, openmpi, or no)"
            ${EXIT_CMD} 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\""
        echo "       The MPI mode is required (choose mpich,  openmpi, or no)"
        ${EXIT_CMD} 1
      fi
      shift 2
      ;;
    -rc | --rebuild_cp2k)
      REBUILD_CP2K="yes"
      shift 1
      ;;
    -t | --test)
      RUN_TEST="yes"
      if (($# > 1)); then
        TESTOPTS="${2}"
      else
        echo "ERROR: No argument found for flag \"${1}\""
        echo "       A string argument with the TESTOPTS (even an empty one \"\") is required"
        ${EXIT_CMD} 1
      fi
      shift 2
      ;;
    -ue | --use_externals)
      USE_EXTERNALS="yes"
      shift 1
      ;;
    -v | --verbose)
      VERBOSE=1
      INSTALL_MESSAGE="LAZY"
      VERBOSE_FLAG="--verbose"
      VERBOSE_MAKEFILE="ON"
      shift 1
      ;;
    --)
      shift 1
      break
      ;;
    -*)
      echo "ERROR: Unknown option \"${1}\" specified"
      ${EXIT_CMD} 1
      ;;
    *)
      break
      ;;
  esac
done

# Remove leading zeros from NUM_PROCS
NUM_PROCS=$(awk '{print $1+0}' <<< "${NUM_PROCS}")

# Check if we are working within a docker or podman container
[[ -f /.dockerenv || -f /run/.containerenv ]] && IN_CONTAINER="yes" || IN_CONTAINER="no"

# Assemble CMake feature flag list
CMAKE_FEATURE_FLAGS="${CMAKE_FEATURE_FLAG_ALL} ${CMAKE_FEATURE_FLAG_MPI} ${CMAKE_FEATURE_FLAGS}"

# Clean CMake feature flag list from repeated entries
declare -A seen=()
out=()
for flag in ${CMAKE_FEATURE_FLAGS}; do
  [[ ${seen[${flag}]+_} ]] || {
    seen[${flag}]=1
    out+=("${flag}")
  }
done
CMAKE_FEATURE_FLAGS="$(printf '%s\n' "${out[*]}")"

export BUILD_DEPS BUILD_DEPS_ONLY BUILD_TYPE CMAKE_FEATURE_FLAGS CMAKE_FEATURE_FLAGS_GPU CRAY CUDA_SM_CODE \
  DISABLE_LOCAL_CACHE GCC_VERSION GPU_MODEL HAS_PODMAN IN_CONTAINER INSTALL_MESSAGE MPI_MODE NUM_PROCS \
  REBUILD_CP2K RUN_TEST TESTOPTS VERBOSE VERBOSE_FLAG VERBOSE_MAKEFILE

# Show help if requested
if [[ "${HELP}" == "yes" ]]; then
  echo ""
  echo "Usage: ${SCRIPT_NAME} [-bd | --build_deps]"
  echo "                    [-bd_only | --build_deps_only]"
  echo "                    [-bt | --build_type (Debug | Release | RelWithDebInfo)]"
  echo "                    [-cray]"
  echo "                    [-cv | --cp2k_version (psmp | ssmp | ssmp-static)]"
  echo "                    [-df | --disable | --disable_feature (all | FEATURE | PACKAGE | none)"
  echo "                    [-dlc | --disable_local_cache]"
  echo "                    [-ef | --enable | --enable_feature (all | FEATURE | PACKAGE | none)"
  echo "                    [-gm | -gpu  | --gpu_model (<CUDA SM code> | P100 | V100 | T400 | A100 | H100 | H200 | GH200 | none)]"
  echo "                    [-gv | --gcc_version (10 | 11 | 12 | 13 | 14 | 15 | 16)]"
  echo "                    [-h | --help]"
  echo "                    [-ip | --install_path PATH]"
  echo "                    [-j #PROCESSES]"
  echo "                    [-mpi | --mpi_mode (mpich | no | openmpi)]"
  echo "                    [-rc | --rebuild_cp2k]"
  echo "                    [-t | -test \"TESTOPTS\"]"
  echo "                    [-ue | --use_externals]"
  echo "                    [-v | --verbose]"
  echo ""
  echo "Flags:"
  echo " --build_deps         : Force a rebuild of all CP2K dependencies from scratch (removes the spack folder)"
  echo " --build_deps_only    : Rebuild ONLY the CP2K dependencies from scratch (removes the spack folder)"
  echo " --build_type         : Set preferred CMake build type (default: \"Release\")"
  echo " --cp2k_version       : CP2K version to be built (default: \"psmp\")"
  echo " -cray                : Use Cray specific spack configuration"
  echo " --disable_local_cache: Don't add local Spack cache"
  echo " --enable_feature     : Enable feature or package (default: all)"
  echo " --disable_feature    : Disable feature or package"
  echo " --help               : Print this help information"
  echo " --gcc_version        : Use the specified GCC version (default: automatically decided by spack)"
  echo " --gpu_model          : Select GPU model (default: none)"
  echo " --install_path       : Define the CP2K installation path (default: ./install)"
  echo " -j                   : Number of processes used in parallel"
  echo " --mpi_mode           : Set preferred MPI mode (default: \"mpich\")"
  echo " --rebuild_cp2k       : Rebuild CP2K: removes the build folder (default: no)"
  echo " --test               : Perform a regression test run after a successful build"
  echo " --use_externals      : Use external packages installed on the host system. This results in much"
  echo "                        faster build times, but it can also cause conflicts with outdated packages"
  echo "                        pulled in from the host system, e.g. old python or gcc versions"
  echo " --verbose            : Write verbose output"
  echo ""
  echo "Hints:"
  echo " - Remove the folder ${CP2K_ROOT}/build to (re)build CP2K from scratch"
  echo "   (see also --rebuild_cp2k flag)"
  echo " - Remove the folder ${CP2K_ROOT}/spack to (re)build CP2K and all its dependencies from scratch"
  echo "   (see also --build_deps flag)"
  echo " - The folder ${CP2K_ROOT}/install is updated after each successful run"
  echo ""
  echo "Packages: all | ace | cosma | deepmd | dftd4 | dlaf | elpa | fftw3 | greenx | hdf5 | libint2 |"
  echo "          libsmeagol | libtorch | libvdwxc | libxsmm | mimic | openpmd | pexsi | plumed | sirius |"
  echo "          spfft | spglib | spla | tblite | trexio | vori "
  echo ""
  echo "Features: cray_pm_accel_energy | cusolver_mp | dbm_gpu | elpa_gpu | grid_gpu | pw_gpu |"
  echo "          spla_gemm_offloading | unified_memory"
  echo ""
  ${EXIT_CMD}
fi

echo ""
echo "BUILD_DEPS          = ${BUILD_DEPS}"
echo "BUILD_DEPS_ONLY     = ${BUILD_DEPS_ONLY}"
echo "CMAKE_BUILD_TYPE    = ${BUILD_TYPE}"
echo "CP2K_VERSION        = ${CP2K_VERSION}"
echo "CRAY                = ${CRAY}"
echo "DISABLE_LOCAL_CACHE = ${DISABLE_LOCAL_CACHE}"
echo "GCC_VERSION         = ${GCC_VERSION}"
if ((CUDA_SM_CODE > 0)); then
  echo "GPU                 = ${GPU_MODEL} (CUDA SM code: ${CUDA_SM_CODE})"
else
  echo "GPU                 = ${GPU_MODEL}"
fi
echo "INSTALL_PREFIX      = ${INSTALL_PREFIX}"
echo "INSTALL_MESSAGE     = ${INSTALL_MESSAGE}"
echo "IN_CONTAINER        = ${IN_CONTAINER}"
echo "MPI_MODE            = ${MPI_MODE}"
echo "NUM_PROCS           = ${NUM_PROCS} (processes)"
echo "Physical cores      = $(lscpu -p=Core,Socket | grep -v '#' | sort -u | wc -l) (host view)"
echo "REBUILD_CP2K        = ${REBUILD_CP2K}"
echo "RUN_TEST            = ${RUN_TEST}"
if [[ "${RUN_TEST}" == "yes" ]]; then
  echo "TESTOPTS            = \"${TESTOPTS}\""
fi
echo "USE_EXTERNALS       = ${USE_EXTERNALS}"
echo "VERBOSE_FLAG        = ${VERBOSE_FLAG}"
echo "VERBOSE_MAKEFILE    = ${VERBOSE_MAKEFILE}"
echo "VERBOSE             = ${VERBOSE}"
if (($# > 0)); then
  echo "Remaining args   =" "$@" "(not used)"
fi
echo ""
echo "LD_LIBRARY_PATH     = ${LD_LIBRARY_PATH:-}"
echo ""
echo "PATH                = ${PATH:-}"
echo ""
echo "CMAKE_FEATURE_FLAGS = ${CMAKE_FEATURE_FLAGS}"
echo ""

((VERBOSE > 0)) && echo "SED_PATTERN_LIST    = ${SED_PATTERN_LIST}"

# Check if a valid number of processes is requested
if ((NUM_PROCS < 1)); then
  echo "ERROR: The requested number of processes should be larger than 0, found \"${NUM_PROCS}\""
  ${EXIT_CMD} 1
elif ((MAX_PROCS > 0)) && ((NUM_PROCS > MAX_PROCS)); then
  echo "WARNING: The requested number of processes (${NUM_PROCS}) is larger than the detected number of CPU cores (${MAX_PROCS})"
fi

# Check if a valid CMake build type is selected
case "${BUILD_TYPE}" in
  Debug | Release | RelWithDebInfo)
    true
    ;;
  *)
    echo "ERROR: Invalid CMake build type \"${BUILD_TYPE}\" selected"
    ${EXIT_CMD} 1
    ;;
esac

# Check if a valid MPI type is selected
case "${MPI_MODE}" in
  mpich | openmpi)
    if [[ "${CP2K_VERSION}" == "ssmp"* ]]; then
      echo "ERROR: MPI type \"${MPI_MODE}\" specified for building a serial CP2K binary"
      ${EXIT_CMD} 1
    fi
    ;;
  no)
    if [[ "${CP2K_VERSION}" == "psmp" ]]; then
      echo "ERROR: MPI type \"${MPI_MODE}\" specified for building an MPI-parallel CP2K binary"
      ${EXIT_CMD} 1
    fi
    ;;
  *)
    echo "ERROR: Invalid MPI type \"${MPI_MODE}\" selected"
    echo "       Choose from (mpich | no | openmpi)"
    ${EXIT_CMD} 1
    ;;
esac

# Check if CP2K_VERSION and the selected features are compatible
case "${CP2K_VERSION}" in
  ssmp | ssmp-static)
    for package in cosma dlaf elpa libfabric libsmeagol mimic openpmd pexsi plumed sirius spla; do
      if [[ "${CMAKE_FEATURE_FLAGS}" == *" -DCP2K_USE_${package^^}=ON"* ]]; then
        echo -e "ERROR: The feature ${package^^} is not available for building serial CP2K binaries (${CP2K_VERSION})\n"
        ${EXIT_CMD} 1
      fi
    done
    # Further exclusions are needed for statically linked serial CP2K binaries
    if [[ "${CP2K_VERSION}" == "ssmp-static" ]]; then
      for package in ace deepmd greenx hdf5 libtorch trexio; do
        if [[ "${CMAKE_FEATURE_FLAGS}" == *" -DCP2K_USE_${package^^}=ON"* ]]; then
          echo -e "ERROR: The feature ${package^^} is not available for building statically linked serial CP2K binaries (${CP2K_VERSION})\n"
          ${EXIT_CMD} 1
        fi
      done
    fi
    ;;
esac

# Perform CUDA GPU related settings
if ((CUDA_SM_CODE > 0)); then
  if command -v nvcc &> /dev/null; then
    CUDA_VERSION=$(nvcc -V | sed -n 's/.*release \([0-9.]*\).*/\1/p')
    NVCC_VERSION=$(awk -v x="${CUDA_VERSION}" 'BEGIN{print 10*x}')
  else
    echo -e "\nERROR: No CUDA toolkit installation found (nvcc compiler not found)\n"
    ${EXIT_CMD} 1
  fi
  # Check if the selected CUDA SM code is valid when the nvidia-smi command is available
  if command -v nvidia-smi &> /dev/null; then
    echo -e "NVIDIA driver installation found:\n"
    nvidia-smi
    HOST_CUDA_SM_CODE=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | tail -n 1 | awk '{print 10*$1}')
    if ((CUDA_SM_CODE > HOST_CUDA_SM_CODE)); then
      echo ""
      echo "ERROR: The requested CUDA SM code (${CUDA_SM_CODE}) is larger than the maximum"
      echo "       CUDA SM code (${HOST_CUDA_SM_CODE}) supported by the host system"
      echo ""
      ${EXIT_CMD} 1
    fi
  else
    echo "INFO: No NVIDIA driver installation found (nvidia-smi not found)"
  fi
  # Check GPU support by the CUDA SDK
  if ((NVCC_VERSION < 128)) && ((CUDA_SM_CODE > 90)); then
    echo ""
    echo "ERROR: The CUDA SDK version ${NVCC_VERSION} does not support GPUs with a"
    echo "       CUDA SM code larger than 90 (found ${CUDA_SM_CODE})"
    echo ""
    ${EXIT_CMD} 1
  fi
  if ((NVCC_VERSION > 129)) && ((CUDA_SM_CODE < 75)); then
    echo ""
    echo "ERROR: The CUDA SDK version ${NVCC_VERSION} does not support GPUs with a"
    echo "       CUDA SM code smaller than 75 (found ${CUDA_SM_CODE})"
    echo ""
    ${EXIT_CMD} 1
  fi
  CMAKE_CUDA_FLAGS="-DCP2K_USE_ACCEL=CUDA"
  CMAKE_CUDA_FLAGS+=" -DCP2K_WITH_GPU=${GPU_MODEL}"
  CMAKE_CUDA_FLAGS+=" -DCMAKE_CUDA_ARCHITECTURES=${CUDA_SM_CODE}"
  CMAKE_CUDA_FLAGS+=" ${CMAKE_FEATURE_FLAGS_GPU}"
  out=()
  for flag in ${CMAKE_CUDA_FLAGS}; do
    [[ ${seen[${flag}]+_} ]] || {
      seen[${flag}]=1
      out+=("${flag}")
    }
  done
  CMAKE_CUDA_FLAGS="$(printf '%s\n' "${out[*]}")"
  echo -e "\nCMAKE_CUDA_FLAGS    = ${CMAKE_CUDA_FLAGS}"
  [[ -n "${CUDA_VERSION:-}" ]] && echo -e "\nCUDA_VERSION        = ${CUDA_VERSION} (nvcc compiler)"
  [[ -n "${CUDA_HOME:-}" ]] && echo -e "\nCUDA_HOME           = ${CUDA_HOME}"
  echo ""
else
  CMAKE_CUDA_FLAGS="-DCP2K_USE_ACCEL=OFF"
fi
export CMAKE_CUDA_FLAGS

### Build CP2K dependencies with Spack if needed or requested ###

# Spack version
export SPACK_VERSION="${SPACK_VERSION:-1.1.1}"
export SPACK_BUILD_PATH="${CP2K_ROOT}/spack"
export SPACK_ROOT="${SPACK_BUILD_PATH}/spack"

# Define the CP2K spack configuration file
export CP2K_CONFIG_FILE="${SPACK_BUILD_PATH}/cp2k_deps_${CP2K_VERSION}.yaml"

# If requested, remove the spack folder for (re)building all CP2K dependencies
if [[ "${BUILD_DEPS}" == "always" ]]; then
  for folder in ${SPACK_BUILD_PATH} ${CP2K_ROOT}/build; do
    if [[ -d "${folder}" ]]; then
      echo "Removing folder \"${folder}\""
      rm -rf "${folder}"
    fi
  done
fi

# If requested, remove the build folder for (re)building CP2K
if [[ "${REBUILD_CP2K}" == "yes" ]]; then
  for folder in ${CP2K_ROOT}/build; do
    if [[ -d "${folder}" ]]; then
      echo "Removing folder \"${folder}\""
      rm -rf "${folder}"
    fi
  done
fi

if [[ ! -d "${SPACK_BUILD_PATH}" ]]; then

  # Create a new local spack folder
  mkdir -p "${SPACK_BUILD_PATH}"
  cd "${SPACK_BUILD_PATH}" || ${EXIT_CMD} 1

  # Reset the spack environment
  unset SPACK_ENV

  # Install Spack
  if [[ ! -d "${SPACK_ROOT}" ]]; then
    echo "Installing Spack ${SPACK_VERSION}"
    wget -q "https://github.com/spack/spack/archive/v${SPACK_VERSION}.tar.gz"
    tar -xzf "v${SPACK_VERSION}.tar.gz" && rm -f "v${SPACK_VERSION}.tar.gz"
    mv -f "${SPACK_BUILD_PATH}/spack-${SPACK_VERSION}" "${SPACK_ROOT}"
  fi
  export PATH="${SPACK_ROOT}/bin:${PATH}"

  # Isolate user configuration for spack
  export SPACK_DISABLE_LOCAL_CONFIG=true
  export SPACK_USER_CONFIG_PATH="${SPACK_BUILD_PATH}"
  export SPACK_USER_CACHE_PATH="${SPACK_BUILD_PATH}/cache"
  mkdir -p "${SPACK_USER_CACHE_PATH}"

  if [[ "${DISABLE_LOCAL_CACHE}" == "no" ]]; then
    # Create and activate a virtual environment (venv) for Python packages
    if command -v python3 -m venv --help &> /dev/null; then
      echo "Installing virtual environment for Python packages"
      if ! python3 -m venv "${SPACK_BUILD_PATH}/venv"; then
        echo "ERROR: The creation of a virtual environment (venv) for Python packages failed"
        ${EXIT_CMD} 1
      fi
      export PATH="${SPACK_BUILD_PATH}/venv/bin:${PATH}"
    else
      echo "ERROR: python3 -m venv was not found"
      ${EXIT_CMD} 1
    fi
    # Upgrade pip and install boto3
    if command -v python3 -m pip --version &> /dev/null; then
      python3 -m pip install "${VERBOSE_FLAG}" --upgrade pip
      echo "Installing boto3 module"
      python3 -m pip install "${VERBOSE_FLAG}" boto3==1.38.11 google-cloud-storage==3.1.0
    else
      echo "ERROR: python3 -m pip was not found"
      ${EXIT_CMD} 1
    fi
  fi

  # Initialize Spack shell hooks
  # shellcheck source=/dev/null
  source "${SPACK_ROOT}/share/spack/setup-env.sh"

  if [[ "${IN_CONTAINER}" == "no" ]] && [[ "${DISABLE_LOCAL_CACHE}" == "no" ]]; then
    # The package podman is required for using a MinIO cache
    if command -v podman &> /dev/null; then
      # Check podman version
      ((VERBOSE > 0)) && podman version
      PODMAN_VERSION=$(podman version --format '{{.Client.Version}}' | awk -F'[.-]' '{print $1}')
      if ((PODMAN_VERSION < 4)); then
        echo "WARNING: Outdated podman version $(podman version --format '{{.Client.Version}}') found"
      fi
      HAS_PODMAN="yes"
    else
      echo "INFO: podman was not found"
      echo "      Install the package podman to take advantage of a local spack cache"
      echo "      This accelerates a rebuild of the CP2K dependencies, significantly"
    fi
  fi

  # Start Spack cache if we are not within a container and have podman available
  if [[ "${IN_CONTAINER}" == "no" ]] && [[ "${DISABLE_LOCAL_CACHE}" == "no" ]]; then
    if [[ "${HAS_PODMAN}" == "yes" ]]; then
      if ! "${CP2K_ROOT}"/tools/docker/spack_cache_start.sh; then
        echo "ERROR: Could not start (new) spack cache"
        echo ""
        echo "An error message starting with \"Error: initial journal cursor: ...\" indicates that the"
        echo "journald logging does not work. Try to switch podman (3.x) to file-based logs by creating"
        echo "the file \"~/.config/containers/containers.conf\" with the following two lines:"
        echo "[containers]"
        echo "log_driver = \"k8s-file\""
        ${EXIT_CMD} 1
      fi
    fi
  fi

  # Add local spack cache if possible
  if [[ "${DISABLE_LOCAL_CACHE}" == "no" ]]; then
    export SPACK_CACHE=${SPACK_CACHE:-"s3://spack-cache --s3-endpoint-url=http://localhost:9000"}
    echo "SPACK_CACHE = \"${SPACK_CACHE}\""
    spack mirror list
    if ! spack mirror list | grep -q "local-cache"; then
      echo "Setting up local spack cache"
      if ((VERBOSE > 0)); then
        "${CP2K_ROOT}"/tools/docker/scripts/setup_spack_cache.sh
      else
        "${CP2K_ROOT}"/tools/docker/scripts/setup_spack_cache.sh &> /dev/null
      fi
    fi
  else
    echo "INFO: A local Spack cache is NOT used"
  fi

  # Prepare the CP2K spack configuration file
  sed -E \
    -e "s|root: /opt/spack|root: ${SPACK_ROOT}/opt/spack|" \
    -e "/\"build_type=/s|build_type=[^\"]*|build_type=${BUILD_TYPE}|" \
    "${CP2K_ROOT}/tools/spack/cp2k_deps_${CP2K_VERSION}.yaml" > "${CP2K_CONFIG_FILE}"

  # Apply selected MPI type if needed
  # MPICH is selected by default for psmp and no change needed for ssmp
  if [[ "${MPI_MODE}" == "openmpi" ]]; then
    sed -E \
      -e '/\s*-\s+"mpich@/ s/^ /#/' \
      -e '/\s*#\s*-\s+"openmpi@/ s/#/ /' \
      -e '/\s*-\s+mpich/ s/mpich$/openmpi/' \
      -i "${CP2K_CONFIG_FILE}"
  fi

  # Activate CUDA in the spack configuration file if requested
  if ((CUDA_SM_CODE > 0)); then
    sed -E \
      -e "0,/~cuda/s//+cuda cuda_arch=${CUDA_SM_CODE}/" \
      -e 's/"~cuda\s+~gpu_direct"/"\+cuda \+gpu_direct"/' \
      -e '/\s*#\s*-\s+"fabrics=efa,ucx"/ s/#/ /' \
      -i "${CP2K_CONFIG_FILE}"
    # Building libfabric with CUDA causes problems
    # sed -E -e 's/"~cuda\s+~gdrcopy"/"\+cuda \+gdrcopy"/' -i "${CP2K_CONFIG_FILE}"
    sed -E -e 's/"~cuda\s+~gdrcopy"/"\~cuda"/' -i "${CP2K_CONFIG_FILE}"
    if [[ -n "${CUDA_VERSION:-}" ]]; then
      # Set CUDA SM code
      sed -E -e "s/spec:\s+cuda@[.0-9]*/spec: cuda@${CUDA_VERSION}/" -i "${CP2K_CONFIG_FILE}"
    fi
    if [[ -n "${CUDA_HOME:-}" ]]; then
      sed -E -e "s|prefix: /usr/local/cuda|prefix: ${CUDA_HOME}|" -i "${CP2K_CONFIG_FILE}"
    fi
  else
    sed -E -e 's/"~cuda\s+~gdrcopy"/"\~cuda"/' -i "${CP2K_CONFIG_FILE}"
  fi

  # Apply Cray specific adaptation of the spack configuration if requested (CSCS)
  if [[ "${CRAY}" == "yes" ]]; then
    sed -E \
      -e '/\s*#\s*-\s+"netmod=ofi"/ s/#/ /' \
      -e 's/~xpmem/+xpmem/' \
      -e 's/"libfabric@[.0-9]*"/"libfabric@1.22.0"/' \
      -i "${CP2K_CONFIG_FILE}"
  fi

  # Apply feature selection to spack configuration file
  if [[ -n "${SED_PATTERN_LIST}" ]]; then
    eval sed -E "${SED_PATTERN_LIST}" -i "${CP2K_CONFIG_FILE}"
  fi

  # Find all compilers
  echo "Searching for GCC compilers"
  if ! spack compiler find &> /dev/null; then
    echo "ERROR: The compiler detection of spack failed"
    ${EXIT_CMD} 1
  fi
  spack compiler list

  # Retrieve the newest compiler version found by spack
  GCC_VERSION_NEWEST="$(spack compilers | awk '/gcc/ {print $2}' | sort -V | tail -n 1)"
  echo "The newest GCC compiler version found by spack is ${GCC_VERSION_NEWEST}"
  GCC_VERSION_NEWEST="$(echo "${GCC_VERSION_NEWEST}" | sed -E -e 's/.*@([0-9]+).*/\1/' | cut -d. -f1)"

  # Check if the newest compiler version found on the host system is new enough
  GCC_VERSION_MINIMUM="10"
  if ((GCC_VERSION_NEWEST < GCC_VERSION_MINIMUM)); then
    echo "INFO: The newest GCC compiler version ${GCC_VERSION_NEWEST} is too old,"
    echo "      because at least version ${GCC_VERSION_MINIMUM} is required"
    GCC_VERSION="14"
    echo "INFO: The lastest stable GCC version ${GCC_VERSION} will be built and used by spack"
  fi

  # Add a specific GCC version to the spack configuration if requested
  if [[ "${GCC_VERSION}" == "auto" ]]; then
    echo "Spack will automatically select the GCC version which is not necessarily the newest one found"
  else
    sed -E -e "s/gcc@10:/gcc@${GCC_VERSION}/" -i "${CP2K_CONFIG_FILE}"
  fi

  # Disable PEXSI because of an issue with SuperLU using recent GCC versions
  if ((CUDA_SM_CODE > 0)) || ((GCC_VERSION_NEWEST > 14)); then
    sed -E -e '/\s*-\s+"pexsi@/ s/^ /#/' -i "${CP2K_CONFIG_FILE}"
    echo "INFO: PEXSI has been disabled because CUDA or GCC 15 is used"
  fi

  # Create CP2K environment if needed
  if spack env list | grep -q "${CP2K_ENV}"; then
    if [[ -n "${SPACK_ENV}" ]]; then
      echo "The Spack environment \"${CP2K_ENV}\" exists already"
    fi
  else
    cat "${CP2K_CONFIG_FILE}"
    echo "The Spack environment \"${CP2K_ENV}\" does NOT exist"
    if spack env create "${CP2K_ENV}" "${CP2K_CONFIG_FILE}" &> /dev/null; then
      echo "The Spack environment \"${CP2K_ENV}\" has been created"
    else
      echo "ERROR: The creation of the Spack environment \"${CP2K_ENV}\" failed"
      ${EXIT_CMD} 1
    fi
  fi

  # Activate CP2K environment if needed
  if spack env activate "${CP2K_ENV}" &> /dev/null; then
    echo "The Spack environment \"${CP2K_ENV}\" has been activated"
  else
    echo "ERROR: The activation of the Spack environment \"${CP2K_ENV}\" failed"
    ${EXIT_CMD} 1
  fi

  # Update the repo builtin
  if ! spack repo update; then
    echo "ERROR: The update of the repo builtin failed"
    ${EXIT_CMD} 1
  fi

  # Add the local CP2K development Spack repository when missing
  export CP2K_REPO="cp2k_dev"
  if ! spack repo list | grep -q "${CP2K_REPO}"; then
    spack repo add --scope "env:${CP2K_ENV}" "${CP2K_ROOT}/tools/spack/spack_repo/${CP2K_REPO}"
  else
    echo "The CP2K development repo ${CP2K_REPO} is already present ... skipping add"
  fi

  spack -e ${CP2K_ENV} repo list

  # Find all external packages
  if [[ "${USE_EXTERNALS}" == "yes" ]]; then
    if ! spack -C "${SPACK_USER_CONFIG_PATH}" external find; then
      echo -e "\nERROR: The detection of externals by spack failed"
      ${EXIT_CMD} 1
    fi
    # Concretize CP2K dependencies
    if ! spack -e "${CP2K_ENV}" concretize --fresh; then
      echo -e "\nERROR: The spack concretize for environment \"${CP2K_ENV}\" failed"
      echo ""
      echo "HINT: The (-ue | --use_externals) flags can cause conflicts with outdated"
      echo "      packages on the host system, e.g. old python or gcc versions"
      echo ""
      ${EXIT_CMD} 1
    fi
  else
    # If the installation of a specific GCC version for building CP2K is requested,
    # then an installation of the same GCC version on the host system has to be removed
    # from the configuration for a static build of CP2K
    if [[ "${CP2K_VERSION}" == *"-static" ]]; then
      if spack compilers | grep -q "gcc@${GCC_VERSION}"; then
        echo "INFO: A static build of CP2K with GCC ${GCC_VERSION} is requested without using"
        echo "      any externals which requires to build that GCC version with spack"
        echo "      while removing a gcc@${GCC_VERSION} installed on the host from the spack"
        echo "      configuration at the same time"
        if ! spack compiler remove "gcc@${GCC_VERSION}"; then
          echo -e "\nERROR: Removing the gcc@${GCC_VERSION} compiler failed"
          ${EXIT_CMD} 1
        fi
      fi
    fi
    if ! spack -e "${CP2K_ENV}" --no-user-config --no-system-config concretize --fresh; then
      echo -e "\nERROR: The spack concretize for environment \"${CP2K_ENV}\" failed"
      ${EXIT_CMD} 1
    fi
  fi

  ((VERBOSE > 0)) && spack find -c

  # Create spack makefile for all dependencies
  if ! spack -e "${CP2K_ENV}" env depfile -o spack_makefile; then
    echo "ERROR: The creation of the spack makefile failed"
    ${EXIT_CMD} 1
  fi

  # Install CP2K dependencies via Spack
  if ! make -j"${NUM_PROCS}" --file=spack_makefile SPACK_COLOR=never --output-sync=recurse; then
    echo "ERROR: Building the CP2K dependencies with spack failed"
    if [[ "${USE_EXTERNALS}" == "yes" ]]; then
      echo "HINT:  Try to re-run the build without the (-ue | --use_externals) flag which avoids"
      echo "       errors or conflicts caused by externals from the host system"
    fi
    ${EXIT_CMD} 1
  fi

  # Find all compilers
  if ! spack compiler find; then
    echo "ERROR: The compiler detection of spack failed"
    ${EXIT_CMD} 1
  fi

  # Fix Libs list in elpa pkg-config file when CUDA is used
  if ((CUDA_SM_CODE > 0)); then
    ELPA_PKG_CONFIG_FILE="$(find -L "${SPACK_ROOT}/opt/spack/view" -name "elpa*.pc")"
    if [[ -f "${ELPA_PKG_CONFIG_FILE}" ]]; then
      sed -E -e 's|Libs: |Libs: -L/usr/local/cuda/lib64 |' -i "${ELPA_PKG_CONFIG_FILE}"
    fi
  fi

  # Return from spack folder after all installations are done
  cd "${CP2K_ROOT}" || ${EXIT_CMD} 1

  # Make a note of the successful build
  touch "${SPACK_BUILD_PATH}/BUILD_DEPENDENCIES_COMPLETED"

  echo -e '\n*** Installation of CP2K dependencies completed ***\n'

else

  # Check if the CP2K dependencies have been built successfully
  if [[ ! -f "${SPACK_BUILD_PATH}/BUILD_DEPENDENCIES_COMPLETED" ]]; then
    echo "ERROR: The last build of the CP2K dependencies was not completed successfully"
    echo "       Re-run the script with the \"--build_dependencies\" or \"-bd\" flag or"
    echo "       remove the folder ${SPACK_BUILD_PATH}"
    ${EXIT_CMD} 1
  fi

  # Initialize Spack shell hooks
  # shellcheck source=/dev/null
  source "${SPACK_ROOT}"/share/spack/setup-env.sh

fi

# Quit after (re)building all CP2K dependencies if requested
if [[ "${BUILD_DEPS_ONLY}" == "yes" ]]; then
  ${EXIT_CMD}
fi

### End of build CP2K dependencies ###

### Build CP2K ###

if spack env list | grep -q "${CP2K_ENV}"; then
  spack env list
else
  echo "ERROR: No Spack environment \"${CP2K_ENV}\" found"
  ${EXIT_CMD} 1
fi

# Activate spack environment
eval "$(spack env activate --sh ${CP2K_ENV})"

spack env status

# CMake configuration step
export CMAKE_BUILD_PATH="${CP2K_ROOT}/build"

# PyTorch's TorchConfig.cmake is buried in the Python site-packages directory
Torch_DIR="$(dirname "$(find "${SPACK_ROOT}" ! -type l -name TorchConfig.cmake | tail -n 1)")"
export Torch_DIR

# Check if PEXSI was built
if [[ "${MPI_MODE}" != "no" ]]; then
  CP2K_USE_PEXSI="$(grep -Eq '\s*#\s*-\s+"pexsi@' spack/cp2k_deps_psmp.yaml && echo OFF || echo ON)"
else
  CP2K_USE_PEXSI="ON"
fi
export CP2K_USE_PEXSI

if [[ ! -d "${CMAKE_BUILD_PATH}" ]]; then
  mkdir -p "${CMAKE_BUILD_PATH}"
  case "${CP2K_VERSION}" in
    "psmp")
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" \
        -GNinja \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH=ON \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        ${CMAKE_FEATURE_FLAGS} \
        -DCP2K_USE_PEXSI="${CP2K_USE_PEXSI}" \
        ${CMAKE_CUDA_FLAGS} \
        -Werror=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    "ssmp")
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" \
        -GNinja \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH=ON \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        ${CMAKE_FEATURE_FLAGS} \
        ${CMAKE_CUDA_FLAGS} \
        -Werror=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    "ssmp-static")
      # Find some static libraries in advance
      LIBOPENBLAS=$(find -L "${SPACK_ROOT}"/opt/spack/view -name libopenblas.a)
      LIBM="$(find /usr -name libm.a 2> /dev/null)"
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" \
        -GNinja \
        -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_EXE_LINKER_FLAGS="-static" \
        -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH=ON \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        -DCP2K_BLAS_LINK_LIBRARIES="${LIBOPENBLAS};${LIBM}" \
        -DCP2K_LAPACK_LINK_LIBRARIES="${LIBOPENBLAS};${LIBM}" \
        ${CMAKE_FEATURE_FLAGS} \
        -Werror=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      # It is almost impossible to avoid that shared libraries are pulled in
      # when the CP2K dependencies are built and thus just change all shared
      # to static libraries hoping that these are also available
      sed -E -e "s/\.so/.a/g" -i build/build.ninja
      ;;
    *)
      echo "ERROR: Unknown CP2K version ${CP2K_VERSION} specified"
      ${EXIT_CMD} 1
      ;;
  esac
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The CMake configuration step failed with the error code ${EXIT_CODE}"
    echo "       You can try to remove the build folder with 'rm -rf build' and re-run"
    echo "       or even start the whole CP2K installation from scratch with the \"-bd\" flag"
    ${EXIT_CMD} "${EXIT_CODE}"
  fi
fi

# CMake build step
echo -e '\n*** Compiling CP2K ***\n'
cmake --build "${CMAKE_BUILD_PATH}" --parallel "${NUM_PROCS}" -- "${VERBOSE_FLAG}" |& tee "${CMAKE_BUILD_PATH}"/ninja.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  echo "ERROR: The CMake build step failed with the error code ${EXIT_CODE}"
  ${EXIT_CMD} "${EXIT_CODE}"
fi

# CMake install step
echo -e '\n*** Installing CP2K ***\n'
cmake --install "${CMAKE_BUILD_PATH}" |& tee "${CMAKE_BUILD_PATH}"/install.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  echo -e "\nERROR: The CMake installation step failed with the error code ${EXIT_CODE}"
  ${EXIT_CMD} "${EXIT_CODE}"
fi

# Collect and compress all log files when building within a container
if [[ "${IN_CONTAINER}" == "yes" ]]; then
  if ! cat "${CMAKE_BUILD_PATH}"/cmake.log \
    "${CMAKE_BUILD_PATH}"/ninja.log \
    "${CMAKE_BUILD_PATH}"/install.log |
    gzip > "${CP2K_ROOT}"/install/build_cp2k.log.gz; then
    echo -e "\nERROR: The compressed log file generation failed"
    ${EXIT_CMD} 1
  fi
fi

# Cut extension from the CP2K version (e.g. ssmp-static -> ssmp)
export VERSION=${CP2K_VERSION%%-*}

# Add CP2K lib folder to the LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${INSTALL_PREFIX}/lib:${LD_LIBRARY_PATH}"

# Assemble search paths for libraries
search_paths="${SPACK_ROOT}/opt/spack/view ${INSTALL_PREFIX}/lib"

# Retrieve paths to all libraries needed by the CP2k binary
for library in $(readelf -d "${INSTALL_PREFIX}/bin/cp2k.${VERSION}" "${INSTALL_PREFIX}"/lib/libcp2k.* | awk '/NEEDED/{print $5}' | tr -d '[]' | sort | uniq); do
  # Search for missing library path
  if ! ldconfig -p | grep -qE "/${library}$"; then
    # shellcheck disable=SC2086
    library_with_path=$(find -L ${search_paths} -name "${library}" | tail -n 1)
    if [[ -n "${library_with_path}" ]]; then
      library_path="$(dirname "${library_with_path}")"
      case ":${LD_LIBRARY_PATH}:" in
        *":${library_path}:"*)
          ((VERBOSE > 1)) && echo "The library path ${library_path} is already present ... skipping"
          ;;
        *)
          LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${library_path}"
          ((VERBOSE > 0)) && echo "The library path ${library_path} was appended to LD_LIBRARY_PATH"
          ;;
      esac
    fi
  fi
done
export LD_LIBRARY_PATH
echo -e "\nLD_LIBRARY_PATH = \"${LD_LIBRARY_PATH}\""

# Create an ld configuration file for CP2K because LD_LIBRARY_PATH is fragile
if [[ "${IN_CONTAINER}" == "yes" ]]; then
  CP2K_LD_CONF_FILE="/etc/ld.so.conf.d/cp2k.conf"
else
  CP2K_LD_CONF_FILE="${INSTALL_PREFIX}/bin/cp2k.conf"
fi

# Either read LD_LIBRARY_PATH from an existing configuration file or write a new one
if [[ -f "${CP2K_LD_CONF_FILE}" ]]; then
  # Load LD_LIBRARY_PATH from an existing configuration file
  LD_LIBRARY_PATH="$(grep -E '^/' "${CP2K_LD_CONF_FILE}" | paste -sd ':' -)"
else
  # Write LD_LIBRARY_PATH to a new configuration file
  echo "${LD_LIBRARY_PATH}" | tr ':' '\n' | grep '^/' > "${CP2K_LD_CONF_FILE}"
fi

export LD_LIBRARY_PATH

# Create links to CP2K binaries
cd "${INSTALL_PREFIX}"/bin || ${EXIT_CMD} 1
for binary in *."${VERSION}"; do
  if ! ln -sf "${binary}" "${binary%%."${VERSION}"}"; then
    echo "ERROR: The creation of a symbolic link for the binary ${binary} failed"
    ${EXIT_CMD}
  fi
done
ln -sf cp2k."${VERSION}" cp2k."${VERSION/smp/opt}"
ln -sf cp2k."${VERSION}" cp2k_shell
cd "${CP2K_ROOT}" || ${EXIT_CMD} 1

# Allow to run as root with OpenMPI
if [[ "${MPI_MODE}" == "openmpi" ]]; then
  OMPI_VARS="export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_plm_rsh_agent=/bin/false"
else
  OMPI_VARS=""
fi
export OMPI_VARS

# MPICH and OpenMPI have different flags for exporting environment variables to MPI ranks
if [[ "${MPI_MODE}" == "openmpi" ]]; then
  ENV_VAR_FLAG="-x"
else
  ENV_VAR_FLAG="-genv"
fi
export ENV_VAR_FLAG

# Assemble flags for running the regression tests
TESTOPTS="--cp2kdatadir ${INSTALL_PREFIX}/share/cp2k/data  --maxtasks ${NUM_PROCS} --workbasedir ${INSTALL_PREFIX}/regtesting ${TESTOPTS}"
export TESTOPTS

# Create launch script for environment setup
LAUNCH_SCRIPT="${INSTALL_PREFIX}"/bin/launch
export LAUNCH_SCRIPT
cat << *** > "${LAUNCH_SCRIPT}"
#!/bin/bash
ulimit -c 0 -s unlimited
export PATH=${INSTALL_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-2}
export OMP_STACKSIZE=256M
${OMPI_VARS}
exec "\$@"
***
chmod 750 "${LAUNCH_SCRIPT}"

# Create shortcut for launching the regression tests
cat << *** > "${INSTALL_PREFIX}"/bin/run_tests
ldd -- ${INSTALL_PREFIX}/bin/cp2k.${VERSION} 2>&1 | grep -E 'not ' | sort | uniq
${CP2K_ROOT}/tests/do_regtest.py ${TESTOPTS} \$* ${INSTALL_PREFIX}/bin ${VERSION}
***
chmod 750 "${INSTALL_PREFIX}"/bin/run_tests

# Optionally, launch test run
if [[ "${RUN_TEST}" == "yes" ]]; then
  echo -e "\n*** Launching regression test run using the script ${INSTALL_PREFIX}/bin/run_tests\n"
  ${LAUNCH_SCRIPT} run_tests
  EXIT_CODE=$?
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The regression test run failed with the error code ${EXIT_CODE}"
    ${EXIT_CMD} "${EXIT_CODE}"
  fi
else
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    if ((CUDA_SM_CODE > 0)); then
      DEVICE_FLAG=" --device nvidia.com/gpu=all"
    else
      DEVICE_FLAG=""
    fi
    echo ""
    echo "*** A regression test run can be launched with"
    echo "    podman run -it${DEVICE_FLAG} --rm <IMAGE ID> run_tests"
    echo ""
    if [[ "${VERSION}" == "ssmp" ]]; then
      echo "*** A CP2K run using 8 OpenMP threads (default) can be launched with"
      echo "    podman run -it${DEVICE_FLAG} --rm <IMAGE ID> cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    else
      echo "*** An MPI-parallel CP2K run using 2 OpenMP threads for each of the 4 MPI ranks can be launched with"
      echo "    podman run -it${DEVICE_FLAG} --rm <IMAGE ID> mpiexec -n 4 ${ENV_VAR_FLAG} OMP_NUM_THREADS=2 cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    fi
    echo ""
  else
    echo ""
    echo "*** A regression test run can be launched with"
    echo "    ${LAUNCH_SCRIPT} run_tests"
    echo ""
    if [[ "${VERSION}" == "ssmp" ]]; then
      echo "*** A CP2K run using 8 OpenMP threads (default) can be launched with"
      echo "    ${LAUNCH_SCRIPT} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
      echo ""
      echo "*** A CP2K run using only 4 OpenMP threads can be launched with"
      echo "    export OMP_NUM_THREADS=4; ${LAUNCH_SCRIPT} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    else
      echo "*** An MPI-parallel CP2K run using 2 OpenMP threads for each of the 4 MPI ranks can be launched with"
      echo "    export OMP_NUM_THREADS=2; ${LAUNCH_SCRIPT} mpiexec -n 4 cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    fi
    echo ""
  fi
fi
