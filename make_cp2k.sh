#!/usr/bin/env bash

# Purpose: Build CP2K using Spack and CMake locally within the folder CP2K_ROOT
#          which defaults to the current working directory. This should be a
#          "cp2k/" folder containing the CP2K source tree which can be downloaded with
#
#          > git clone https://github.com/cp2k/cp2k.git cp2k
#
#          and after
#
#          > cd cp2k
#
#          this script can be run using the default options with
#
#          > ./make_cp2k.sh
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
#          It is recommended to store already compiled packages in a cache (default).
#          This will accelerate the (re)build of the CP2K dependencies with Spack
#          significantly (see -uc flag).
#
#          After the CP2K dependencies are built with Spack, CP2K itself is built
#          and installed using CMake in the subfolders cp2k/build and cp2k/install,
#          respectively.
#
#          Subsequent runs of the script will use the CMake configuration in the
#          subfolder cp2k/build. A rebuild of CP2K from scratch can be enforced
#          by removing or renaming that subfolder (see also -rc flag).
#
#          The folders cp2k/spack and cp2k/build are the default locations within
#          the build path (see -bp flag) and cp2k/install is the default installation
#          path (see -ip flag).
#
#          A CP2K regression run can be launched automatically by adding the flag
#          -t "" (or --test ""). This flag expects a string with the TESTOPTS, e.g.
#
#          > ./make_cp2k.sh -t "--maxtasks 8 --restrictdir QS/regtest-gpw-1"
#
#          Alternatively, the script cp2k/install/bin/run_tests can be launched
#          after a successful CP2K build. Usage examples are printed at the end.
#
#          The work is performed in four steps:
#
#          (1) make_cp2k.sh: Parse the flags, check the settings, remove folders
#              if requested (-bd, -rc), and write the resolved configuration to
#              cp2k/install/config/options.conf or only print it (--print_config)
#          (2) tools/make_cp2k/build_deps.sh: Build the CP2K dependencies with
#              Spack if needed and write cp2k/install/config/deps.conf
#          (3) tools/make_cp2k/build_cp2k.sh: CMake configuration, build, and
#              installation of CP2K and write cp2k/install/config/build.conf
#          (4) tools/make_cp2k/{check_conventions,run_regtests,test_coverage,
#              test_performance,build_gromacs,build_ase}.sh:
#              Optional checks, tests, and builds as requested by the flags
#              (test_performance.sh also creates the run_benchmarks script for psmp)
#
#          The step scripts are executed (not sourced) and read their settings
#          only from the configuration files passed as arguments. Thus each step
#          can also be re-run on its own, e.g. a rebuild of CP2K with
#
#          > tools/make_cp2k/build_cp2k.sh install/config/{options,deps}.conf

# Authors: Matthias Krack (MK)

# Version: 3.0

# The script has to be executed and must not be sourced, because it would
# otherwise modify the calling shell (shell options, variables, exit)
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
  echo "ERROR: $(basename "${BASH_SOURCE[0]}") must be executed, not sourced"
  return 1
fi

# Treat unset variables as an error and let a pipeline fail if any of its commands fails
set -uo pipefail

### Step 1: Parse the flags, check the settings, and write the resolved configuration ###

# Retrieve script name
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"

# Load the helper functions shared with the step scripts
MAKE_CP2K_PATH="$(dirname "$(realpath "${BASH_SOURCE[0]}")")/tools/make_cp2k"
if [[ ! -f "${MAKE_CP2K_PATH}/common.sh" ]]; then
  echo "ERROR: The file ${MAKE_CP2K_PATH}/common.sh was not found"
  exit 1
fi
# shellcheck source=/dev/null
source "${MAKE_CP2K_PATH}/common.sh"

# Pre-scan the arguments for the flag --print_config, because it has to be known
# before any output is written: all regular output is then redirected to stderr
# and only the resolved configuration is written to stdout
PRINT_CONFIG="no"
for arg in "$@"; do
  case "${arg}" in
    -pc | --print_config | --print-config)
      PRINT_CONFIG="yes"
      ;;
  esac
done
if [[ "${PRINT_CONFIG}" == "yes" ]]; then
  # Keep stdout as file descriptor 3 for the configuration output
  exec 3>&1 1>&2
fi

# Check platform
if [[ "$(uname -s)" == "Darwin" ]]; then
  echo "ERROR: A build for Apple macOS (Darwin, Homebrew) is not supported yet,"
  echo "       but CP2K can be built within a container using podman."
  echo "       Check the Dockerfiles in the folder cp2k/tools/docker like"
  echo "       Dockerfile.test_spack_openmpi-psmp and the usage description in their header"
  exit 1
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
    exit 1
  fi
done

# Disable python bytecode generation to suppress __pycache__ folder creation
export PYTHONDONTWRITEBYTECODE=1

# Abort git transfers slower than 1 KB/s for 60 s, since git has no timeout,
# so that a stalled clone (e.g. py-torch with submodules) fails and is retried
export GIT_HTTP_LOW_SPEED_LIMIT=1000
export GIT_HTTP_LOW_SPEED_TIME=60

# Print OS info if available
if [[ -f /etc/os-release ]]; then
  echo ""
  cat /etc/os-release
fi

# Check bash version
if ((BASH_VERSINFO < 4)); then
  echo ""
  echo "ERROR: The employed bash version ${BASH_VERSION} is too old (from 2004)"
  echo "       Install a newer bash version"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    echo "HINT:  Use the bash command provided by Homebrew with Apple"
    echo "       /opt/homebrew/bin should precede /bin/bash in the PATH"
  fi
  exit 1
elif ((BASH_VERSINFO < 5)); then
  echo -e "\nWARNING: The employed bash version ${BASH_VERSION} is quite old (from 2009)"
else
  echo -e "\nINFO: Using bash version ${BASH_VERSION}"
fi

# Check if the python3 version is new enough for spack
if ! python3 -c "import sys; sys.exit(not(sys.version_info >= (3, 10)))"; then
  echo ""
  echo "ERROR: Python version is NOT >= 3.10 (needed for Spack)"
  echo "       Found only $(python3 -V)"
  exit 1
else
  echo -e "\nINFO: Found $(python3 -V)"
fi

# Convert a CMake build type to its canonical spelling (case-insensitive input),
# unknown build types are returned unchanged
canonical_build_type() {
  case "${1,,}" in
    debug) echo "Debug" ;;
    release) echo "Release" ;;
    relwithdebinfo) echo "RelWithDebInfo" ;;
    *) echo "${1}" ;;
  esac
}

# Remove repeated CMake flags (-DNAME=VALUE) from a list. The last definition
# of each flag is kept at the position of its first appearance.
unique_cmake_flags() {
  local flag name
  local -A last=()
  local -a order=() out=()
  for flag in ${1}; do
    name="${flag%%=*}"
    [[ -n "${last[${name}]+_}" ]] || order+=("${name}")
    last[${name}]="${flag}"
  done
  for name in "${order[@]+"${order[@]}"}"; do
    out+=("${last[${name}]}")
  done
  echo "${out[*]+"${out[*]}"}"
}

# Default values
ASE_VERSION=""
BENCHMARK_PROFILE=""
BUILD_DEPS="if_needed"
BUILD_DEPS_ONLY="no"
BUILD_SHARED_LIBS="${BUILD_SHARED_LIBS:-ON}"
CHECK_CONVENTIONS="no"
CP2K_BUILD_TYPE="$(canonical_build_type "${CP2K_BUILD_TYPE:-Release}")"
DEPS_BUILD_TYPE="$(canonical_build_type "${DEPS_BUILD_TYPE:-Release}")"
CMAKE_FEATURE_FLAG_ALL="-DCP2K_USE_EVERYTHING=ON" # all features are activated by default
CMAKE_FEATURE_FLAGS="-DCP2K_BLAS_VENDOR=OpenBLAS" # LAPACK/BLAS from OpenBLAS by default
CMAKE_FEATURE_FLAG_MPI="-DCP2K_USE_MPI=ON"        # MPI is switched on by default
CMAKE_FEATURE_FLAGS_GPU="-DCP2K_USE_SPLA_GEMM_OFFLOADING=ON"
CMAKE_PRESET="native-gnu-x86_64"
CRAY="no"
CUDA_SM_CODE=0
GCC_VERSION="auto"
GPU_MODEL="none"
GROMACS_VERSION=""
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
NUM_PACKAGES=2
NVCC_VERSION=0
REBUILD_CP2K="no"
RUN_BENCHMARK="no"
RUN_TEST="no"
SED_PATTERN_LIST=""
TESTOPTS=""
TEST_ASE="no"
TEST_COVERAGE="no"
TEST_GROMACS="no"
USE_CACHE="folder"
USE_CUSOLVER_MP=""
USE_EXTERNALS="no"
USE_OPENCL="no"
VERBOSE=0
VERBOSE_FLAG="--quiet"
VERBOSE_MAKEFILE="OFF"
VERBOSE_SPACK=""

export CP2K_ENV="cp2k_env"
export CP2K_ROOT=${CP2K_ROOT:-${PWD}}
export CP2K_VERSION="${CP2K_VERSION:-psmp}"

# Retrieve CP2K revision if folder is a git repository
if CP2K_REVISION=$(git -C "${CP2K_ROOT}" rev-parse --short HEAD 2> /dev/null); then
  # Git succeeded, CP2K_REVISION is already set
  :
else
  CP2K_REVISION="unknown"
fi

export BUILD_PATH="${BUILD_PATH:-${CP2K_ROOT}}"
export INSTALL_PREFIX="${INSTALL_PREFIX:-${CP2K_ROOT}/install}"

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -ase)
      if (($# > 1)); then
        ASE_VERSION="${2}"
      else
        echo "ERROR: No ASE version (branch or tag name) found for flag \"${1}\""
        exit 1
      fi
      shift 2
      ;;
    -bd | --build_deps | --build_dependencies)
      BUILD_DEPS="always"
      shift 1
      ;;
    -bd_only | --build_deps_only | --build_dependencies_only)
      BUILD_DEPS="always"
      BUILD_DEPS_ONLY="yes"
      shift 1
      ;;
    -bp | --build_path)
      if (($# > 1)); then
        BUILD_PATH="$(realpath "${2}")"
      else
        echo "ERROR: No build path argument found for flag \"${1}\""
        exit 1
      fi
      shift 2
      ;;
    -bsl | --build_static_libcp2k)
      BUILD_SHARED_LIBS="OFF"
      shift 1
      ;;
    -bt | --build_type)
      if (($# > 1)); then
        CP2K_BUILD_TYPE="$(canonical_build_type "${2}")"
        case "${CP2K_BUILD_TYPE}" in
          Debug)
            if [[ "${CP2K_VERSION}" =~ smp$ ]]; then
              CP2K_VERSION="${CP2K_VERSION/smp/dbg}"
              DEPS_BUILD_TYPE="RelWithDebInfo"
            fi
            ;;
          Release | RelWithDebInfo)
            if [[ "${CP2K_VERSION}" =~ dbg$ ]]; then
              CP2K_VERSION="${CP2K_VERSION/dbg/smp}"
              DEPS_BUILD_TYPE="Release"
            fi
            ;;
          *)
            echo "ERROR: Invalid CMake build type \"${2}\" specified (choose Debug, Release, or RelWithDebInfo)"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No CMake build type argument found for flag \"${1}\" (e.g. Release)"
        exit 1
      fi
      shift 2
      ;;
    -cc | --check_conventions)
      CHECK_CONVENTIONS="yes"
      shift 1
      ;;
    -cray)
      CRAY="yes"
      shift 1
      ;;
    -cv | --cp2k_version)
      if (($# > 1)); then
        case "${2,,}" in
          pdbg | psmp | sdbg | ssmp | ssmp-static)
            CP2K_VERSION="${2,,}"
            # Set build type
            case "${CP2K_VERSION}" in
              pdbg | sdbg)
                CP2K_BUILD_TYPE="Debug"
                DEPS_BUILD_TYPE="RelWithDebInfo"
                ;;
              psmp | ssmp | ssmp-static)
                CP2K_BUILD_TYPE="Release"
                DEPS_BUILD_TYPE="Release"
                ;;
            esac
            # Disable MPI for a serial CP2K binary
            case "${CP2K_VERSION}" in
              sdbg | ssmp | ssmp-static)
                MPI_MODE="no"
                CMAKE_FEATURE_FLAG_MPI="-DCP2K_USE_MPI=OFF"
                # LibGint requires MPI
                CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_LIBGINT=OFF"
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"libgint@/ s/^ /#/'"
                ;;
            esac
            # Apply specific package selection for statically linked binaries
            case "${CP2K_VERSION}" in
              ssmp-static)
                CMAKE_FEATURE_FLAG_ALL="-DCP2K_USE_EVERYTHING=ON"
                for package in libfci libint2 libxc libxs spglib vori tblite; do
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_${package^^}=ON"
                done
                for package in ace deepmd gauxc greenx hdf5 libtorch pexsi trexio; do
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_${package^^}=OFF"
                done
                ;;
            esac
            ;;
          *)
            echo "ERROR: Invalid CP2K version \"${2}\" specified (choose pdbg, psmp, sdbg, ssmp, or ssmp-static)"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\" (choose pdbg, psmp, sdbg, ssmp, or ssmp-static)"
        exit 1
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
            # Enable or disable all features
            CMAKE_FEATURE_FLAG_ALL="-DCP2K_USE_EVERYTHING=${ON_OFF}"
            for package in adios2 cosma deepmdkit dla-future dla-future-fortran elpa \
              gauxc greenx hdf5 libfabric libfci libint2 libvdwxc libsmeagol libvori \
              libxc libxs libxsmm mimic-mcl openpmd-api pace pexsi plumed py-torch sirius \
              spfft spglib spla tblite trexio; do
              SED_PATTERN_LIST+=" -e '/\s*-\s+\"${package}@/ ${SUBST}"
            done
            # dbcsr must use blas as fallback when libxs/libxsmm is disabled
            if [[ "${ON_OFF}" == "OFF" ]]; then
              SED_PATTERN_LIST+=" -e '/\s*-\s+\"smm=libxs\"/ s/libxs/blas/'"
            fi
            ;;
          ace | cosma | deepmd | dftd4 | dlaf | elpa | gauxc | greenx | hdf5 | libfci | libgint | \
            libint2 | libsmeagol | libtorch | libxc | libxs | mimic | openpmd | pexsi | plumed | \
            spglib | tblite | trexio | vori)
            CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_${2^^}=${ON_OFF}"
            # Translate package selection to sed pattern
            case "${2,,}" in
              ace)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"p${2,,}@/ ${SUBST}"
                ;;
              cosma | elpa | greenx | hdf5 | libfci | libgint | libsmeagol | libxc | pexsi | plumed | \
                spglib | trexio)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}@/ ${SUBST}"
                ;;
              deepmd)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}kit@/ ${SUBST}"
                if [[ "${ON_OFF}" == "ON" ]]; then
                  # DeePMD-kit requires libtorch
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_LIBTORCH=${ON_OFF}"
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"py-torch@/ ${SUBST}"
                fi
                ;;
              dftd4)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}@/ ${SUBST}"
                if [[ "${ON_OFF}" == "ON" ]]; then
                  echo "INFO: tblite is disabled to avoid conflicts when dftd4 is used"
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_TBLITE=OFF"
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"tblite@/ s/^ /#/'"
                fi
                ;;
              dlaf)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"dla-future.*@/ ${SUBST}"
                if [[ "${ON_OFF}" == "ON" ]]; then
                  SED_PATTERN_LIST+=" -e 's/\~dlaf/\+dlaf/'"
                else
                  SED_PATTERN_LIST+=" -e 's/\+dlaf/\~dlaf/'"
                fi
                ;;
              gauxc)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}@/ ${SUBST}"
                if [[ "${ON_OFF}" == "ON" ]]; then
                  # GauXC requires libtorch
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_LIBTORCH=${ON_OFF}"
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"py-torch@/ ${SUBST}"
                fi
                ;;
              libint2)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"libint@/ ${SUBST}"
                ;;
              libtorch)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"py-torch@/ ${SUBST}"
                if [[ "${ON_OFF}" == "OFF" ]]; then
                  # DeePMD-kit and GauXC require libtorch
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_DEEPMD=${ON_OFF}"
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"deepmdkit@/ ${SUBST}"
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_GAUXC=${ON_OFF}"
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"gauxc@/ ${SUBST}"
                fi
                ;;
              libxs)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"libxs@/ ${SUBST}"
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"libxsmm@/ ${SUBST}"
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
              tblite)
                SED_PATTERN_LIST+=" -e '/\s*-\s+\"${2,,}@/ ${SUBST}"
                if [[ "${ON_OFF}" == "ON" ]]; then
                  echo "INFO: dftd4 is disabled to avoid conflicts when tblite is used"
                  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_DFTD4=OFF"
                  SED_PATTERN_LIST+=" -e '/\s*-\s+\"dftd4@/ s/^ /#/'"
                fi
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
          cusolver_mp)
            USE_CUSOLVER_MP="${ON_OFF}"
            CMAKE_FEATURE_FLAGS_GPU+=" -DCP2K_USE_${2^^}=${ON_OFF}"
            SED_PATTERN_LIST+=" -e '/\s*-\s+\"cusolvermp@/ ${SUBST}"
            ;;
          cray_pm_accel_energy | spla_gemm_offloading | unified_memory)
            CMAKE_FEATURE_FLAGS_GPU+=" -DCP2K_USE_${2^^}=${ON_OFF}"
            ;;
          dbm_gpu | elpa_gpu | grid_gpu | libxc_gpu | fft_gpu)
            CMAKE_FEATURE_FLAGS_GPU+=" -DCP2K_ENABLE_${2^^}=${ON_OFF}"
            ;;
          none)
            # do nothing
            ;;
          *)
            echo "ERROR: Unknown CP2K feature \"${2}\" specified"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No feature found for flag \"${1}\""
        exit 1
      fi
      shift 2
      ;;
    -gm | -gpu | --gpu_model)
      if (($# > 1)); then
        case "${2^^}" in
          P100 | V100 | T400 | A100 | A40 | H100 | H200 | GH200 | B200)
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
              B200)
                CUDA_SM_CODE=100
                ;;
            esac
            ;;
          60 | 70 | 75 | 80 | 86 | 87 | 89 | 90 | 100 | 120 | 121)
            CUDA_SM_CODE=${2}
            ;;
          NONE)
            GPU_MODEL="${2,,}"
            ;;
          *)
            echo -e "\nERROR: Unknown GPU model \"${2}\" specified (choose <CUDA SM code>, P100, V100, T400, A100, A40, H100, H200, GH200, B200 or none)\n"
            exit 1
            ;;
        esac
      else
        echo -e "\nERROR: No argument found for flag \"${1}\" (choose <CUDA SM code>, P100, V100, T400, A100, A40, H100, H200, GH200, B200 or none)\n"
        exit 1
      fi
      shift 2
      ;;
    -gromacs)
      if (($# > 1)); then
        GROMACS_VERSION="${2}"
      else
        echo "ERROR: No GROMACS version (branch or tag name) found for flag \"${1}\" (e.g. v2026.3)"
        exit 1
      fi
      shift 2
      ;;
    -gv | --gcc_version)
      if (($# > 1)); then
        case "${2}" in
          1[0-6])
            GCC_VERSION="${2}"
            ;;
          *)
            echo "ERROR: Invalid GCC version \"${2}\" specified (choose from 10 to 16)"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\" (choose a GCC version)"
        exit 1
      fi
      shift 2
      ;;
    -h | --help)
      HELP="yes"
      shift 1
      ;;
    -ip | --install_path | --install_prefix)
      if (($# > 1)); then
        INSTALL_PREFIX="$(realpath "${2}")"
      else
        echo "ERROR: No install path argument found for flag \"${1}\""
        exit 1
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
            exit 1
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
            exit 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\""
        echo "       The MPI mode is required (choose mpich,  openmpi, or no)"
        exit 1
      fi
      shift 2
      ;;
    -np | --num_packages)
      if (($# > 1)); then
        case "${2}" in
          -*)
            shift 1
            ;;
          [0-9]*)
            NUM_PACKAGES="${2}"
            shift 2
            ;;
          *)
            echo "ERROR: The -np flag can only be followed by an integer number, found \"${2}\""
            exit 1
            ;;
        esac
      else
        shift 1
      fi
      ;;
    -opencl)
      USE_OPENCL="yes"
      shift 1
      ;;
    -pc | --print_config | --print-config)
      # Already handled by the pre-scan of the arguments
      shift 1
      ;;
    -preset)
      if (($# > 1)); then
        case "${2}" in
          asan | native-gnu-x86_64 | native-gnu-arm64 | native-intel | none)
            CMAKE_PRESET="${2}"
            ;;
          *)
            echo "ERROR: Invalid preset \"${2}\" specified"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No CMake preset found for flag \"${1}\""
        exit 1
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
        exit 1
      fi
      shift 2
      ;;
    -ta | --test_ase)
      TEST_ASE="yes"
      shift 1
      ;;
    -tc | --test_coverage)
      TEST_COVERAGE="yes"
      shift 1
      ;;
    -tg | --test_gromacs)
      TEST_GROMACS="yes"
      shift 1
      ;;
    -tp | --test_performance)
      RUN_BENCHMARK="yes"
      if (($# > 1)); then
        case "${2,,}" in
          cuda_* | default | openmp)
            BENCHMARK_PROFILE="${2,,}"
            ;;
          *)
            echo "ERROR: Invalid benchmark profile \"${2}\" specified (choose e.g. openmp)"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No benchmark profile found for flag \"${1}\""
        exit 1
      fi
      shift 2
      ;;
    -uc | --use_cache)
      if (($# > 1)); then
        case "${2,,}" in
          folder | minio | no | none)
            USE_CACHE="${2,,}"
            ;;
          *)
            echo "ERROR: Invalid cache mode \"${2}\" specified (choose folder, minio or no)"
            exit 1
            ;;
        esac
      else
        echo "ERROR: No argument found for flag \"${1}\" (choose folder, minio or no)"
        exit 1
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
      VERBOSE_SPACK="--verbose"
      shift 1
      ;;
    --)
      shift 1
      break
      ;;
    -*)
      echo "ERROR: Unknown option \"${1}\" specified"
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

# Remove leading zeros from NUM_PACKAGES and NUM_PROCS
NUM_PACKAGES=$(awk '{print $1+0}' <<< "${NUM_PACKAGES}")
NUM_PROCS=$(awk '{print $1+0}' <<< "${NUM_PROCS}")

# Check if we are working within a docker or podman container
[[ -f /.dockerenv || -f /run/.containerenv ]] && IN_CONTAINER="yes" || IN_CONTAINER="no"

# Assemble CMake feature flag list
CMAKE_FEATURE_FLAGS="${CMAKE_FEATURE_FLAG_ALL} ${CMAKE_FEATURE_FLAG_MPI} ${CMAKE_FEATURE_FLAGS}"

# DLA-Future does not work with MPICH yet
case "${MPI_MODE}" in
  mpich | openmpi)
    if [[ "${CMAKE_FEATURE_FLAGS}" == *"-DCP2K_USE_EVERYTHING=ON"* ]] ||
      [[ "${CMAKE_FEATURE_FLAGS}" == *"-DCP2K_USE_DLAF=ON"* ]]; then
      echo -e "\nINFO: DLA-Future does not work with ${MPI_MODE^^} yet and is disabled"
      CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_DLAF=OFF"
      SED_PATTERN_LIST+=" -e '/\s*-\s+\"dla-future.*@/ s/^ /#/'"
      SED_PATTERN_LIST+=" -e 's/\+dlaf/\~dlaf/'"
    fi
    ;;
esac

# Set default GROMACS version if GROMACS/CP2K testing is requested
if [[ ${TEST_GROMACS} == "yes" ]]; then
  GROMACS_VERSION=${GROMACS_VERSION:-v2026.3}
fi

# Set default ASE version if ASE/CP2K testing is requested
if [[ ${TEST_ASE} == "yes" ]]; then
  ASE_VERSION=${ASE_VERSION:-master}
fi

# Perform setup for an AddressSanitizer (ASAN) build
if [[ "${CMAKE_PRESET}" == "asan" ]]; then
  if [[ "${CHECK_CONVENTIONS}" == "yes" ]] || [[ "${TEST_COVERAGE}" == "yes" ]]; then
    echo -e "\nERROR: The CMake preset \"asan\" cannot be combined with convention checking or coverage testing\n"
    exit 1
  fi
  if [[ "${CP2K_VERSION}" == "ssmp-static" ]]; then
    echo -e "\nERROR: The CMake preset \"asan\" cannot be used with CP2K_VERSION \"ssmp-static\"\n"
    exit 1
  fi
  # The ASAN flags of the preset are only defined for the build type RelWithDebInfo
  if [[ "${CP2K_BUILD_TYPE}" != "RelWithDebInfo" ]]; then
    echo ""
    echo "WARNING: The CMake preset \"asan\" implies CP2K_BUILD_TYPE \"RelWithDebInfo\" but found \"${CP2K_BUILD_TYPE}\""
    echo "         CP2K_BUILD_TYPE is set to \"RelWithDebInfo\""
    echo ""
    CP2K_BUILD_TYPE="RelWithDebInfo"
    if [[ "${CP2K_VERSION}" =~ dbg$ ]]; then
      CP2K_VERSION="${CP2K_VERSION/dbg/smp}"
      DEPS_BUILD_TYPE="Release"
    fi
  fi
  # TODO Re-enable GREENX. It currently leads to a heap-buffer-overflow
  # in `greenx_refine_pade()` at greenx_interface.F:80.
  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_GREENX=OFF"
fi

# Perform setup for coding conventions check
if [[ "${CHECK_CONVENTIONS}" == "yes" ]]; then
  if [[ "${CP2K_VERSION}" != "psmp" ]]; then
    echo ""
    echo "WARNING: Convention checking implies CP2K_VERSION \"psmp\" but found \"${CP2K_VERSION}\""
    echo "         CP2K_VERSION is set to \"psmp\""
    echo ""
    CP2K_VERSION="psmp"
  fi
  CP2K_BUILD_TYPE="RelWithDebInfo"
  CMAKE_PRESET="conventions"
  Fortran_COMPILER_LAUNCHER="${CP2K_ROOT}/tools/conventions/redirect_gfortran_output.py"
else
  Fortran_COMPILER_LAUNCHER=""
fi

# Test coverage and generate a coverage report
if [[ "${TEST_COVERAGE}" == "yes" ]]; then
  if [[ "${CP2K_VERSION}" != "psmp" ]]; then
    echo ""
    echo "WARNING: Coverage testing implies CP2K_VERSION \"psmp\" but found \"${CP2K_VERSION}\""
    echo "         CP2K_VERSION is set to \"psmp\""
    echo ""
    CP2K_VERSION="psmp"
  fi
  if ! command -v lcov &> /dev/null; then
    echo -e "\nERROR: The package lcov is mandatory for a coverage test"
    echo -e "       Install the missing package and re-run the script\n"
    exit 1
  fi
  CP2K_BUILD_TYPE="RelWithDebInfo"
  CMAKE_PRESET="coverage"
  RUN_TEST="yes"
  TESTOPTS+=" --ompthreads=1 --keepalive"
fi

# Clean the CMake feature flag list from repeated entries
CMAKE_FEATURE_FLAGS="$(unique_cmake_flags "${CMAKE_FEATURE_FLAGS}")"

export ASE_VERSION BENCHMARK_PROFILE BUILD_DEPS BUILD_DEPS_ONLY BUILD_SHARED_LIBS CHECK_CONVENTIONS CMAKE_FEATURE_FLAGS \
  CMAKE_FEATURE_FLAGS_GPU CP2K_BUILD_TYPE CP2K_REVISION CRAY CUDA_SM_CODE DEPS_BUILD_TYPE Fortran_COMPILER_LAUNCHER \
  GCC_VERSION GPU_MODEL GROMACS_VERSION IN_CONTAINER INSTALL_MESSAGE MPI_MODE NUM_PACKAGES NUM_PROCS \
  REBUILD_CP2K RUN_BENCHMARK RUN_TEST TEST_COVERAGE TEST_ASE TEST_GROMACS TESTOPTS USE_CACHE USE_OPENCL VERBOSE \
  VERBOSE_FLAG VERBOSE_MAKEFILE VERBOSE_SPACK

# Show help if requested
if [[ "${HELP}" == "yes" ]]; then
  echo ""
  echo "Usage: ${SCRIPT_NAME} [-bd | --build_deps]"
  echo "                    [-ase ASE_VERSION]"
  echo "                    [-bd_only | --build_deps_only]"
  echo "                    [-bp | --build_path PATH]"
  echo "                    [-bsl | --build_static_libcp2k]"
  echo "                    [-bt | --build_type (Debug | Release | RelWithDebInfo)]"
  echo "                    [-cc | --check_conventions]"
  echo "                    [-cray]"
  echo "                    [-cv | --cp2k_version (pdbg | psmp | sdbg | ssmp | ssmp-static)]"
  echo "                    [-df | --disable | --disable_feature (all | FEATURE | PACKAGE | none)]"
  echo "                    [-ef | --enable | --enable_feature (all | FEATURE | PACKAGE | none)]"
  echo "                    [-gm | -gpu  | --gpu_model (<CUDA SM code> | P100 | V100 | T400 | A100 | H100 | H200 | GH200 | B200 | none)]"
  echo "                    [-gromacs GROMACS_VERSION]"
  echo "                    [-gv | --gcc_version (10 | 11 | 12 | 13 | 14 | 15 | 16)]"
  echo "                    [-h | --help]"
  echo "                    [-ip | --install_path PATH]"
  echo "                    [-j #PROCESSES]"
  echo "                    [-mpi | --mpi_mode (mpich | no | openmpi)]"
  echo "                    [-np | --num_packages #PACKAGES]"
  echo "                    [-opencl]"
  echo "                    [-pc | --print_config]"
  echo "                    [-preset (asan | native-gnu-x86_64 | native-gnu-arm64 | native-intel | none)]"
  echo "                    [-rc | --rebuild_cp2k]"
  echo "                    [-t | --test \"TESTOPTS\"]"
  echo "                    [-ta | --test_ase]"
  echo "                    [-tc | --test_coverage]"
  echo "                    [-tg | --test_gromacs]"
  echo "                    [-tp | --test_performance \"BENCHMARK_PROFILE\"]"
  echo "                    [-uc | --use_cache (folder | minio | no | none)]"
  echo "                    [-ue | --use_externals]"
  echo "                    [-v | --verbose]"
  echo ""
  echo "Flags:"
  echo " -ase                  : Build CP2K with ASE support"
  echo " --build_deps          : Force a rebuild of all CP2K dependencies from scratch (removes the spack folder)"
  echo " --build_deps_only     : Rebuild ONLY the CP2K dependencies from scratch (removes the spack folder)"
  echo " --build_path          : Define the CP2K build path (default: ${CP2K_ROOT})"
  echo " --build_static_libcp2k: Build a static CP2K library libcp2k.a instead of the default shared one libcp2k.so"
  echo " --build_type          : Set preferred CMake build type for CP2K (default: \"Release\")"
  echo " --check_conventions   : Check compliance with CP2K's coding conventions"
  echo " --cp2k_version        : CP2K version to be built (default: \"psmp\")"
  echo " -cray                 : Use Cray specific spack configuration"
  echo " --enable_feature      : Enable feature or package (default: all)"
  echo " --disable_feature     : Disable feature or package"
  echo " -gromacs              : Build GROMACS with CP2K support"
  echo " --help                : Print this help information"
  echo " --gcc_version         : Use the specified GCC version (default: automatically decided by spack)"
  echo " --gpu_model           : Select GPU model (default: none)"
  echo " --install_path        : Define the CP2K installation path (default: ./install)"
  echo " -j                    : Maximum number of processes (CPU cores) used in parallel"
  echo "                         If the variable OMP_NUM_THREADS is set and the -j flag is not supplied, then"
  echo "                         the maximum number of processes is defined by OMP_NUM_THREADS"
  echo " --mpi_mode            : Set preferred MPI mode (default: \"mpich\")"
  echo " --num_packages        : Maximum number of packages built by spack in parallel (default: 2)"
  echo " -opencl               : Enable the use of the Open Computing Language (OpenCL)"
  echo " --print_config        : Print the resolved configuration (options.conf) as KEY=value lines to stdout"
  echo "                         and quit without building anything (all other output goes to stderr)"
  echo " -preset               : Use a CMake configure preset, see \"cmake --list-presets\" (default: native-gnu-x86_64)"
  echo " --rebuild_cp2k        : Rebuild CP2K: removes the build folder (default: no)"
  echo " --test                : Perform a regression test run after a successful build"
  echo " --test_ase            : Build and test CP2K with ASE support"
  echo " --test_coverage       : Analyse the code coverage and generate a coverage report"
  echo " --test_gromacs        : Build and test GROMACS with CP2K support"
  echo " --test_performance    : Perform a benchmark run after a successful build"
  echo " --use_cache           : Use a \"folder\", a \"MinIO\" object storage container (requires podman) or \"no\" cache"
  echo "                         Set the environment variable SPACK_CACHE to specify the folder name, e.g."
  echo "                         SPACK_CACHE=\"file://${CP2K_ROOT}/spack_cache\" (default)"
  echo " --use_externals       : Use external packages installed on the host system. This results in much"
  echo "                         faster build times, but it can also cause conflicts with outdated packages"
  echo "                         pulled in from the host system, e.g. old python or gcc versions"
  echo " --verbose             : Write verbose output"
  echo ""
  echo "Hints:"
  echo " - Remove the folder ${CP2K_ROOT}/build to (re)build CP2K from scratch"
  echo "   (see also --rebuild_cp2k flag)"
  echo " - Remove the folder ${CP2K_ROOT}/spack to (re)build CP2K and all its dependencies from scratch"
  echo "   (see also --build_deps flag)"
  echo " - The folder ${CP2K_ROOT}/install is updated after each successful run"
  echo ""
  echo "Packages: all | ace | cosma | deepmd | dftd4 | dlaf | elpa | gauxc | greenx | hdf5 | libfci |"
  echo "          libgint | libint2 | libsmeagol | libtorch | libvdwxc | libxs | mimic | openpmd | pexsi | plumed |"
  echo "          sirius | spfft | spglib | spla | tblite | trexio | vori "
  echo ""
  echo "Features: cray_pm_accel_energy | cusolver_mp | dbm_gpu | elpa_gpu | grid_gpu | libxc_gpu | fft_gpu |"
  echo "          spla_gemm_offloading | unified_memory"
  echo ""
  exit
fi

echo ""
if [[ -n ${ASE_VERSION} ]]; then
  echo "ASE_VERSION         = ${ASE_VERSION}"
fi
echo "BUILD_DEPS          = ${BUILD_DEPS}"
echo "BUILD_DEPS_ONLY     = ${BUILD_DEPS_ONLY}"
echo "BUILD_PATH          = ${BUILD_PATH}"
echo "BUILD_SHARED_LIBS   = ${BUILD_SHARED_LIBS}"
echo "CHECK_CONVENTIONS   = ${CHECK_CONVENTIONS}"
echo "CMAKE_PRESET        = ${CMAKE_PRESET}"
echo "CP2K_BUILD_TYPE     = ${CP2K_BUILD_TYPE}"
echo "CP2K_REVISION       = ${CP2K_REVISION}"
echo "CP2K_VERSION        = ${CP2K_VERSION}"
echo "CRAY                = ${CRAY}"
echo "DEPS_BUILD_TYPE     = ${DEPS_BUILD_TYPE}"
echo "GCC_VERSION         = ${GCC_VERSION}"
if [[ -n ${GROMACS_VERSION} ]]; then
  echo "GROMACS_VERSION     = ${GROMACS_VERSION}"
fi
if ((CUDA_SM_CODE > 0)); then
  echo "GPU                 = ${GPU_MODEL} (CUDA SM code: ${CUDA_SM_CODE})"
else
  echo "GPU                 = ${GPU_MODEL}"
fi
echo "INSTALL_PREFIX      = ${INSTALL_PREFIX}"
echo "INSTALL_MESSAGE     = ${INSTALL_MESSAGE}"
echo "IN_CONTAINER        = ${IN_CONTAINER}"
echo "MPI_MODE            = ${MPI_MODE}"
echo "NUM_PACKAGES        = ${NUM_PACKAGES} (packages are built by spack concurrently)"
echo "NUM_PROCS           = ${NUM_PROCS} (processes)"
echo "Physical cores      = $(lscpu -p=Core,Socket | grep -v '#' | sort -u | wc -l) (host view)"
echo "REBUILD_CP2K        = ${REBUILD_CP2K}"
echo "RUN_BENCHMARK       = ${RUN_BENCHMARK}"
if [[ "${RUN_BENCHMARK}" == "yes" ]]; then
  case ${CP2K_VERSION} in
    psmp)
      echo "BENCHMARK_PROFILE   = ${BENCHMARK_PROFILE}"
      ;;
    *)
      echo -e "\nERROR: Performance test runs are supported only for CP2K_VERSION \"psmp\", found version \"${CP2K_VERSION}\"\n"
      exit 1
      ;;
  esac
fi
echo "RUN_TEST            = ${RUN_TEST}"
if [[ "${RUN_TEST}" == "yes" ]]; then
  echo "TESTOPTS            = \"${TESTOPTS}\""
fi
echo "TEST_ASE            = ${TEST_ASE}"
echo "TEST_COVERAGE       = ${TEST_COVERAGE}"
echo "TEST_GROMACS        = ${TEST_GROMACS}"
echo "USE_CACHE           = ${USE_CACHE}"
echo "USE_EXTERNALS       = ${USE_EXTERNALS}"
echo "USE_OPENCL          = ${USE_OPENCL}"
echo "VERBOSE_FLAG        = ${VERBOSE_FLAG}"
echo "VERBOSE_MAKEFILE    = ${VERBOSE_MAKEFILE}"
echo "VERBOSE             = ${VERBOSE}"
if (($# > 0)); then
  echo -e "\nRemaining args   =" "$@" "(not used)"
fi
echo ""
echo "LD_LIBRARY_PATH     = ${LD_LIBRARY_PATH:-}"
echo ""
echo "PATH                = ${PATH:-}"
echo ""
echo "CMAKE_FEATURE_FLAGS = ${CMAKE_FEATURE_FLAGS}"
echo ""

((VERBOSE > 0)) && echo -e "SED_PATTERN_LIST    = ${SED_PATTERN_LIST}\n"

# Check if a valid number for the packages to be built by spack in parallel is given
if ((NUM_PACKAGES < 1)); then
  echo -e "\nERROR: The requested number of packages to be built by spack in parallel should be larger than 0, found \"${NUM_PACKAGES}\""
  exit 1
elif ((NUM_PROCS > 0)) && ((NUM_PACKAGES > NUM_PROCS)); then
  echo -e "\nERROR: The requested number of packages to be built in parallel by spack (${NUM_PACKAGES}) is larger than the requested number of processes (${NUM_PROCS})"
  exit 1
fi

# Check if a valid number of processes is requested
if ((NUM_PROCS < 1)); then
  echo "ERROR: The requested number of processes should be larger than 0, found \"${NUM_PROCS}\""
  exit 1
elif ((MAX_PROCS > 0)) && ((NUM_PROCS > MAX_PROCS)); then
  echo "WARNING: The requested number of processes (${NUM_PROCS}) is larger than the detected number of CPU cores (${MAX_PROCS})"
fi

# Check if a valid CMake build type is selected for the dependencies
case "${DEPS_BUILD_TYPE}" in
  Release | RelWithDebInfo)
    true
    ;;
  Debug)
    echo "ERROR: The CMake build type \"${DEPS_BUILD_TYPE}\" is not supported for building the dependencies"
    exit 1
    ;;
  *)
    echo "ERROR: Invalid CMake build type \"${DEPS_BUILD_TYPE}\" selected for building the dependencies"
    exit 1
    ;;
esac

# Check if a valid CMake build type is selected for CP2K
case "${CP2K_BUILD_TYPE}" in
  Debug | Release | RelWithDebInfo)
    true
    ;;
  *)
    echo "ERROR: Invalid CMake build type \"${CP2K_BUILD_TYPE}\" selected for building CP2K"
    exit 1
    ;;
esac
if [[ "${CP2K_VERSION}" == "ssmp-static" ]] && [[ "${CP2K_BUILD_TYPE}" == "Debug" ]]; then
  echo "ERROR: The CMake build type \"${CP2K_BUILD_TYPE}\" is not supported for CP2K_VERSION \"ssmp-static\""
  exit 1
fi

# Check if a valid MPI type is selected
case "${MPI_MODE}" in
  mpich | openmpi)
    if [[ "${CP2K_VERSION}" == "sdbg" || "${CP2K_VERSION}" == "ssmp"* ]]; then
      echo "ERROR: MPI type \"${MPI_MODE}\" specified for building a serial CP2K binary"
      exit 1
    fi
    USE_MPI="ON"
    ;;
  no)
    if [[ "${CP2K_VERSION}" == "pdbg" || "${CP2K_VERSION}" == "psmp" ]]; then
      echo "ERROR: MPI type \"${MPI_MODE}\" specified for building an MPI-parallel CP2K binary"
      exit 1
    fi
    # shellcheck disable=SC2034 # written to options.conf
    USE_MPI="OFF"
    ;;
  *)
    echo "ERROR: Invalid MPI type \"${MPI_MODE}\" selected"
    echo "       Choose from (mpich | no | openmpi)"
    exit 1
    ;;
esac

# cuSOLVERMp requires both CUDA and MPI support.
if [[ "${USE_CUSOLVER_MP}" == "ON" ]]; then
  if [[ "${MPI_MODE}" == "no" ]]; then
    echo -e "ERROR: The feature CUSOLVER_MP is not available for building serial CP2K binaries (${CP2K_VERSION})\n"
    exit 1
  fi
  if ((CUDA_SM_CODE == 0)); then
    echo -e "ERROR: The feature CUSOLVER_MP requires CUDA support (specify --gpu_model)\n"
    exit 1
  fi
fi

# Check if CP2K_VERSION and the selected features are compatible
case "${CP2K_VERSION}" in
  ssmp | ssmp-static)
    for package in cosma dlaf elpa libfabric libgint libsmeagol mimic openpmd pexsi plumed sirius spla; do
      if [[ "${CMAKE_FEATURE_FLAGS}" == *" -DCP2K_USE_${package^^}=ON"* ]]; then
        echo -e "ERROR: The feature ${package^^} is not available for building serial CP2K binaries (${CP2K_VERSION})\n"
        exit 1
      fi
    done
    # Further exclusions are needed for statically linked serial CP2K binaries
    if [[ "${CP2K_VERSION}" == "ssmp-static" ]]; then
      for package in ace deepmd greenx hdf5 libtorch trexio; do
        if [[ "${CMAKE_FEATURE_FLAGS}" == *" -DCP2K_USE_${package^^}=ON"* ]]; then
          echo -e "ERROR: The feature ${package^^} is not available for building statically linked serial CP2K binaries (${CP2K_VERSION})\n"
          exit 1
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
    exit 1
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
      exit 1
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
    exit 1
  fi
  if ((NVCC_VERSION > 129)) && ((CUDA_SM_CODE < 75)); then
    echo ""
    echo "ERROR: The CUDA SDK version ${NVCC_VERSION} does not support GPUs with a"
    echo "       CUDA SM code smaller than 75 (found ${CUDA_SM_CODE})"
    echo ""
    exit 1
  fi
  CMAKE_CUDA_FLAGS="-DCP2K_USE_ACCEL=CUDA"
  CMAKE_CUDA_FLAGS+=" -DCMAKE_CUDA_HOST_COMPILER=$(which g++)"
  CMAKE_CUDA_FLAGS+=" -DCP2K_WITH_GPU=${GPU_MODEL}"
  CMAKE_CUDA_FLAGS+=" -DCMAKE_CUDA_ARCHITECTURES=${CUDA_SM_CODE}"
  CMAKE_CUDA_FLAGS+=" ${CMAKE_FEATURE_FLAGS_GPU}"
  CMAKE_CUDA_FLAGS="$(unique_cmake_flags "${CMAKE_CUDA_FLAGS}")"
  echo -e "\nCMAKE_CUDA_FLAGS    = ${CMAKE_CUDA_FLAGS}"
  [[ -n "${CUDA_VERSION:-}" ]] && echo -e "\nCUDA_VERSION        = ${CUDA_VERSION} (nvcc compiler)"
  [[ -n "${CUDA_HOME:-}" ]] && echo -e "\nCUDA_HOME           = ${CUDA_HOME}"
  echo ""
else
  CMAKE_CUDA_FLAGS="-DCP2K_USE_ACCEL=OFF"
  # CUDA is required for LibGint
  CMAKE_FEATURE_FLAGS+=" -DCP2K_USE_LIBGINT=OFF"
  CMAKE_FEATURE_FLAGS="$(unique_cmake_flags "${CMAKE_FEATURE_FLAGS}")"
  echo -e "\nINFO: LibGint requires CUDA support which is not enabled"
  echo -e "      The CMAKE_FEATURE_FLAGS have been updated to disable LibGint\n"
fi
export CMAKE_CUDA_FLAGS

# Cut extension from the CP2K version (e.g. ssmp-static -> ssmp)
export VERSION=${CP2K_VERSION%%-*}

# Spack version and installation path
export SPACK_VERSION="${SPACK_VERSION:-1.2.2}"
export SPACK_BUILD_PATH="${BUILD_PATH}/spack"
export SPACK_ROOT="${SPACK_BUILD_PATH}/spack"

# Isolate user configuration for spack
export SPACK_DISABLE_LOCAL_CONFIG=true
export SPACK_USER_CONFIG_PATH="${SPACK_BUILD_PATH}"
export SPACK_USER_CACHE_PATH="${SPACK_BUILD_PATH}/cache"

# Define the CP2K spack configuration file
export CP2K_CONFIG_FILE="${SPACK_BUILD_PATH}/cp2k_deps_${CP2K_VERSION:0:1}${CP2K_VERSION:4}.yaml"

# CMake build path
export CMAKE_BUILD_PATH="${BUILD_PATH}/build"

# Folder for the configuration files exchanged with the step scripts
export CONFIG_PATH="${INSTALL_PREFIX}/config"

# Set image tag if available
export IMAGE_TAG=${IMAGE_TAG:-<IMAGE ID>}

# Variables written to the configuration file options.conf (unset variables are skipped)
OPTIONS_CONFIG="ASE_VERSION BENCHMARK_PROFILE BUILD_DEPS BUILD_DEPS_ONLY BUILD_PATH BUILD_SHARED_LIBS
  CHECK_CONVENTIONS CMAKE_BUILD_PATH CMAKE_CUDA_FLAGS CMAKE_FEATURE_FLAGS CMAKE_FEATURE_FLAGS_GPU CMAKE_PRESET
  CONFIG_PATH CP2K_BUILD_TYPE CP2K_CONFIG_FILE CP2K_ENV CP2K_REVISION CP2K_ROOT CP2K_VERSION CRAY CUDA_HOME
  CUDA_SM_CODE CUDA_VERSION DEPS_BUILD_TYPE Fortran_COMPILER_LAUNCHER GCC_VERSION GPU_MODEL GROMACS_VERSION
  IMAGE_TAG IN_CONTAINER INSTALL_MESSAGE INSTALL_PREFIX MPI_MODE NUM_PACKAGES NUM_PROCS NVCC_VERSION
  REBUILD_CP2K RUN_BENCHMARK RUN_TEST SED_PATTERN_LIST SPACK_BUILD_PATH SPACK_CACHE SPACK_DISABLE_LOCAL_CONFIG
  SPACK_ROOT SPACK_USER_CACHE_PATH SPACK_USER_CONFIG_PATH SPACK_VERSION TEST_ASE TEST_COVERAGE TEST_GROMACS
  TESTOPTS USE_CACHE USE_CUSOLVER_MP USE_EXTERNALS USE_MPI USE_OPENCL VERBOSE VERBOSE_FLAG VERBOSE_MAKEFILE
  VERBOSE_SPACK VERSION"

# Print the resolved configuration and quit if requested. This has to happen
# before anything is created or removed.
if [[ "${PRINT_CONFIG}" == "yes" ]]; then
  # shellcheck disable=SC2086
  write_config ${OPTIONS_CONFIG} >&3
  exit 0
fi

# Retrieve and print available resources
free -h
echo ""
ulimit -c 0 -s unlimited
ulimit -a
echo ""

# If requested, remove the spack folder for (re)building all CP2K dependencies
if [[ "${BUILD_DEPS}" == "always" ]]; then
  for folder in ${SPACK_BUILD_PATH} ${CMAKE_BUILD_PATH} ${INSTALL_PREFIX}; do
    if [[ -d "${folder}" ]]; then
      echo "Removing folder \"${folder}\""
      rm -rf "${folder}"
    fi
  done
fi

# If requested, remove the build folder for (re)building CP2K
if [[ "${REBUILD_CP2K}" == "yes" ]]; then
  for folder in ${CMAKE_BUILD_PATH} ${INSTALL_PREFIX}; do
    if [[ -d "${folder}" ]]; then
      echo "Removing folder \"${folder}\""
      rm -rf "${folder}"
    fi
  done
fi

# Write the resolved configuration for the step scripts
mkdir -p "${CONFIG_PATH}"
# shellcheck disable=SC2086
write_config ${OPTIONS_CONFIG} > "${CONFIG_PATH}/options.conf"
echo -e "INFO: The configuration was written to ${CONFIG_PATH}/options.conf\n"

### Step 2: Build CP2K dependencies with Spack if needed or requested ###

"${MAKE_CP2K_PATH}"/build_deps.sh "${CONFIG_PATH}/options.conf" || exit $?

# Quit after (re)building all CP2K dependencies if requested
if [[ "${BUILD_DEPS_ONLY}" == "yes" ]]; then
  exit
fi

### Step 3: CMake configuration, build, and installation of CP2K ###

"${MAKE_CP2K_PATH}"/build_cp2k.sh "${CONFIG_PATH}/options.conf" "${CONFIG_PATH}/deps.conf" || exit $?

### Step 4: Optional checks, tests, and builds based on the CP2K installation ###

# Run a step script with all configuration files and quit if it fails
run_step() {
  "${MAKE_CP2K_PATH}/${1}" "${CONFIG_PATH}/options.conf" "${CONFIG_PATH}/deps.conf" "${CONFIG_PATH}/build.conf" || exit $?
}

# Check coding conventions
if [[ "${CHECK_CONVENTIONS}" == "yes" ]]; then
  run_step check_conventions.sh
fi

# Run regression tests
if [[ "${RUN_TEST}" == "yes" ]]; then
  run_step run_regtests.sh
fi

# Analyse the code coverage of the regression test run
if [[ "${TEST_COVERAGE}" == "yes" ]]; then
  run_step test_coverage.sh
fi

# Create the benchmark script and run the benchmarks if requested (psmp only)
if [[ "${VERSION}" == "psmp" ]]; then
  run_step test_performance.sh
fi

# Build and test GROMACS with CP2K support
if [[ -n "${GROMACS_VERSION}" ]]; then
  run_step build_gromacs.sh
fi

# Build and test ASE with CP2K support
if [[ -n "${ASE_VERSION}" ]]; then
  run_step build_ase.sh
fi
