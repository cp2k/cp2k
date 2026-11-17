# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Default values and parsing of the command line flags
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

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
          dbm_gpu | elpa_gpu | grid_gpu | libxc_gpu | pw_gpu)
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
