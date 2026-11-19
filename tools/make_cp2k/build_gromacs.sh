#!/usr/bin/env bash
#
# build_gromacs.sh: Build and test GROMACS with CP2K support (step 4 of make_cp2k.sh)
#
# Usage: build_gromacs.sh CONFIG_FILE...
#
# The settings are read from the configuration files options.conf, deps.conf,
# and build.conf written by make_cp2k.sh and build_deps.sh.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"
activate_spack_env

# Optionally, build GROMACS/CP2K
if [[ -n "${GROMACS_VERSION}" ]]; then

  # Download GROMACS
  GROMACS_ROOT="${CMAKE_BUILD_PATH}"/gromacs
  [[ -d "${GROMACS_ROOT}" ]] && rm -rf "${GROMACS_ROOT}"
  echo -e "\n*** Downloading GROMACS ${GROMACS_VERSION} ***\n"
  # Spack's OpenSSL libs on LD_LIBRARY_PATH can break the system git-remote-https
  # helper (mismatched libldap), so clone with a clean library path
  if ! LD_LIBRARY_PATH="" retry_download git clone -b "${GROMACS_VERSION}" -c advice.detachedHead=false --depth=1 "${VERBOSE_FLAG}" --single-branch \
    https://gitlab.com/gromacs/gromacs.git "${GROMACS_ROOT}"; then
    echo -e "\nERROR: Downloading GROMACS ${GROMACS_VERSION} failed"
    exit 1
  fi
  cd "${GROMACS_ROOT}" || exit 1
  GROMACS_REVISION="$(git rev-parse --short HEAD)"

  # Retrieve the GROMACS major version from the source tree, since GROMACS_VERSION
  # is a branch or tag name (e.g. v2026.3, release-2026, main)
  GROMACS_VERSION_MAJOR=$(sed -n -E 's/^ *project\( *Gromacs +VERSION +([0-9]+).*/\1/p' CMakeLists.txt)
  if [[ ! ${GROMACS_VERSION_MAJOR} =~ ^[0-9]+$ ]]; then
    echo "ERROR: Could not determine the GROMACS major version from CMakeLists.txt"
    exit 1
  fi
  echo "GROMACS major version: ${GROMACS_VERSION_MAJOR}"

  # CP2K requires MPI_THREAD_MULTIPLE when attaching to an already-initialized MPI environment.
  # GROMACS >= 2027 requests MPI_THREAD_MULTIPLE if the environment variable
  # GMX_USE_MPI_THREAD_MULTIPLE is set, whereas older versions always request
  # MPI_THREAD_FUNNELED and have to be patched
  if ((GROMACS_VERSION_MAJOR >= 2027)); then
    sed -i '/^exec "\$@"$/i export GMX_USE_MPI_THREAD_MULTIPLE=1' "${LAUNCH_SCRIPT}"
  else
    sed -E -e 's/MPI_Init_thread\(argc, argv, MPI_THREAD_FUNNELED,/MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE,/' \
      -i src/gromacs/utility/init.cpp
  fi

  # CMake configuration step for GROMACS
  GROMACS_BUILD_PATH="${GROMACS_ROOT}"/build
  mkdir -p "${GROMACS_BUILD_PATH}"
  echo -e "\n*** Performing CMake configuration for GROMACS ${GROMACS_VERSION} ***\n"
  cmake -S "${GROMACS_ROOT}" -B "${GROMACS_BUILD_PATH}" \
    -DBUILD_SHARED_LIBS="OFF" \
    -DCMAKE_BUILD_TYPE="${CP2K_BUILD_TYPE}" \
    -DCMAKE_CXX_STANDARD="20" \
    -DCMAKE_INSTALL_LIBDIR="lib" \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    -DCP2K_DIR="${INSTALL_PREFIX}"/lib \
    -DGMX_BUILD_OWN_FFTW="ON" \
    -DGMX_CP2K="ON" \
    -DGMX_DOUBLE="OFF" \
    -DGMX_INSTALL_NBLIB_API="OFF" \
    -DGMX_MPI="${USE_MPI}" \
    -DGMXAPI="OFF" \
    -Werror=dev \
    -Wno-error=non-target-directive \
    &> "${GROMACS_BUILD_PATH}/cmake.log"
  EXIT_CODE=$?
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The CMake configuration step for GROMACS failed with the error code ${EXIT_CODE}"
    [[ "${IN_CONTAINER}" == "yes" ]] && mkdir -p /workspace/artifacts/ && cp "${GROMACS_BUILD_PATH}"/*.log /workspace/artifacts/
    tail -n 1000 "${GROMACS_BUILD_PATH}"/cmake.log
    echo -e "\nStatus: FAILED\n"
    exit "${EXIT_CODE}"
  fi

  # CMake build step for GROMACS
  echo -e "\n*** Compiling GROMACS ${GROMACS_VERSION} ***\n"
  cmake --build "${GROMACS_BUILD_PATH}" --parallel "${NUM_PROCS}" --target all qmmm_applied_forces-test &> "${GROMACS_BUILD_PATH}"/make.log
  EXIT_CODE=${PIPESTATUS[0]}
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The CMake build step for GROMACS failed with the error code ${EXIT_CODE}"
    [[ "${IN_CONTAINER}" == "yes" ]] && mkdir -p /workspace/artifacts/ && cp "${GROMACS_BUILD_PATH}"/*.log /workspace/artifacts/
    tail -n 1000 "${GROMACS_BUILD_PATH}"/make.log
    exit "${EXIT_CODE}"
  fi

  # CMake install step for GROMACS
  echo -e "\n*** Installing GROMACS ***\n"
  cmake --install "${GROMACS_BUILD_PATH}" &> "${GROMACS_BUILD_PATH}"/install.log &&
    cp "${GROMACS_BUILD_PATH}"/bin/qmmm_applied_forces-test "${INSTALL_PREFIX}"/bin
  EXIT_CODE=${PIPESTATUS[0]}
  if ((EXIT_CODE != 0)); then
    [[ "${IN_CONTAINER}" == "yes" ]] && mkdir -p /workspace/artifacts/ && cp "${GROMACS_BUILD_PATH}"/*.log /workspace/artifacts/
    echo -e "\nERROR: The CMake installation step for GROMACS failed with the error code ${EXIT_CODE}"
    tail -n 1000 "${GROMACS_BUILD_PATH}"/install.log
    exit "${EXIT_CODE}"
  fi

  # Suppress GROMACS quote and reminder messages
  export GMX_NO_QUOTES=1

  # Print instructions for testing GROMACS/CP2K
  echo ""
  echo "*** The GROMACS/CP2K installation can be tested with"
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    echo "    podman run -it --rm ${IMAGE_TAG} ${LAUNCH_SCRIPT} qmmm_applied_forces-test"
  else
    echo "    ${LAUNCH_SCRIPT} qmmm_applied_forces-test"
  fi
  echo ""

  # Test GROMACS/CP2K installation
  GROMACS_BINARY="gmx"
  [[ ${USE_MPI} == "ON" ]] && GROMACS_BINARY+="_mpi"
  if [[ "${TEST_GROMACS}" == "yes" ]]; then
    if ${LAUNCH_SCRIPT} ${GROMACS_BINARY} --version; then
      echo -e "\n*** Running GROMACS/CP2K QM/MM unit test ***\n"
      if ${LAUNCH_SCRIPT} qmmm_applied_forces-test; then
        echo -e "\nSummary: GROMACS commit ${GROMACS_REVISION} works fine"
        echo -e "Status: OK\n"
      else
        echo -e "\nSummary: Something is wrong with GROMACS commit ${GROMACS_REVISION}"
        echo -e "Status: FAILED\n"
        exit 0
      fi
    else
      echo -e "\nSummary: Something is wrong with GROMACS commit ${GROMACS_REVISION}"
      echo -e "Status: FAILED\n"
      exit 0
    fi
  fi

  # Print usage hint
  echo -e "\n*** See benchmarks/GROMACS/MQAE/README.md for how to run GROMACS/CP2K\n"

fi
