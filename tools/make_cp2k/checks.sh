# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Settings printout, consistency checks, CUDA settings, and available resources
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

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
  out=()
  for flag in ${CMAKE_CUDA_FLAGS}; do
    [[ ${seen[${flag}]+_} ]] || {
      # shellcheck disable=SC2004 # seen is an associative array declared in settings.sh
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

# Retrieve and print available resources
free -h
echo ""
ulimit -c 0 -s unlimited
ulimit -a
echo ""
