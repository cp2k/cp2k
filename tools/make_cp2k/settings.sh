# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Derived settings: CMake feature flags and build modes (ASAN, conventions, coverage)
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

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

# Clean CMake feature flag list from repeated entries and keep only the last definition
declare -A last=()
order=()
# First pass: record the last occurrence of each flag
for flag in ${CMAKE_FEATURE_FLAGS}; do
  name=${flag%=*}
  last[${name}]=${flag}
  order+=("${name}")
done
# Second pass: output only the last definition per flag
# but respecting original order of last appearance
declare -A seen=()
out=()
for name in "${order[@]}"; do
  if [[ -z ${seen[${name}]+_} ]]; then
    seen[${name}]=1
    out+=("${last[${name}]}")
  fi
done
CMAKE_FEATURE_FLAGS="$(printf '%s\n' "${out[*]}")"

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

export ASE_VERSION BENCHMARK_PROFILE BUILD_DEPS BUILD_DEPS_ONLY BUILD_SHARED_LIBS CHECK_CONVENTIONS CMAKE_FEATURE_FLAGS \
  CMAKE_FEATURE_FLAGS_GPU CP2K_BUILD_TYPE CP2K_REVISION CRAY CUDA_SM_CODE DEPS_BUILD_TYPE Fortran_COMPILER_LAUNCHER \
  GCC_VERSION GPU_MODEL GROMACS_VERSION IN_CONTAINER INSTALL_MESSAGE MPI_MODE NUM_PACKAGES NUM_PROCS \
  REBUILD_CP2K RUN_BENCHMARK RUN_TEST TEST_COVERAGE TEST_ASE TEST_GROMACS TESTOPTS USE_CACHE USE_OPENCL VERBOSE \
  VERBOSE_FLAG VERBOSE_MAKEFILE VERBOSE_SPACK
