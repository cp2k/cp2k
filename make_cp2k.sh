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
#          software stacks.
#
#          It is recommended to install podman to take advantage of a spack cache.
#          This will accelerate the build of the CP2K dependencies with Spack
#          significantly.
#
#          After the CP2K dependencies are built with Spack, CP2K itself is built
#          and installed using CMake in the subfolders cp2k/build and cp2k/install.
#
#          Subsequent runs of the script will use the CMake configuration in the
#          subfolder cp2k/build. A rebuild of CP2K from scratch can be enforced
#          by removing or renaming that subfolder.
#
#          A CP2K regression run can be launched automatically by adding the flag
#          -t "" (or --test ""). This flag expects a string with the TESTOPTS, e.g.
#
#          > source ./make_cp2k.sh -t "--restrictdir QS/regtest-gpw-1"
#
#          Alternatively, the script cp2k/install/run_tests can be launched after
#          a successful CP2K build.

# Authors: Matthias Krack (MK)

# Version: 0.2

# History: - Creation (19.12.2025, MK)
#          - Version 0.1: First working version (09.01.2026, MK)
#          - Version 0.2: Add more flags and checks (19.01.2026, MK)

# set -uo pipefail # can be useful for debugging this script

# Trust other scripts being sourced
# shellcheck source=/dev/null

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

# Check platform
if [[ "$(uname -s)" == "Darwin" ]]; then
  echo "ERROR: Apple (Darwin, Homebrew) is not supported yet"
  echo "       Build CP2K within a container, e.g. using Ubuntu"
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

# Default values
BUILD_DEPS="no"
BUILD_TYPE="RelWithDebInfo"
HELP="no"
INSTALL_MESSAGE="NEVER"
if command -v nproc &> /dev/null; then
  NUM_PROCS=${NUM_PROCS:-$(nproc)}
else
  NUM_PROCS=${NUM_PROCS:-8}
fi
RUN_TEST="no"
TESTOPTS=""
VERBOSE=0
VERBOSE_FLAG="--quiet"
VERBOSE_MAKEFILE="OFF"

export CP2K_ENV="cp2k_env"
export CP2K_ROOT=${PWD}
export CP2K_VERSION="${CP2K_VERSION:-psmp}"
export INSTALL_PREFIX="${INSTALL_PREFIX:-${CP2K_ROOT}/install}"

# Parse flags
while [[ $# -gt 0 ]]; do
  case "$1" in
    -bd | --build_deps)
      BUILD_DEPS="yes"
      ;;
    -bt | --build_type)
      BUILD_TYPE="${2}"
      shift 2
      ;;
    -cv | --cp2k_version)
      CP2K_VERSION="${2}"
      shift 2
      ;;
    -h | --help)
      HELP="yes"
      shift 1
      ;;
    -ip | --install_path)
      INSTALL_PREFIX="${2}"
      shift 2
      ;;
    -j | --num_procs)
      NUM_PROCS="${2}"
      shift 2
      ;;
    -t | --test)
      RUN_TEST="yes"
      TESTOPTS="${2}"
      shift 2
      ;;
    -v | --verbose)
      VERBOSE=1
      INSTALL_MESSAGE="LAZY"
      VERBOSE_FLAG="--verbose"
      VERBOSE_MAKEFILE="ON"
      shift 1
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "ERROR: Unknown option ${1} specified"
      ${EXIT_CMD} 1
      ;;
    *)
      break
      ;;
  esac
done

export BUILD_TYPE INSTALL_MESSAGE NUM_PROCS RUN_TEST TESTOPTS VERBOSE VERBOSE_FLAG VERBOSE_MAKEFILE

# Show help if requested
if [[ "${HELP}" == "yes" ]]; then
  echo ""
  echo "Usage: ${SCRIPT_NAME} [-bd|--build_deps] [-bt|--build_type (Debug|Release|RelWithDebInfo|MinSizeRel)]"
  echo "                    [-cv|--cp2k_version (psmp|ssmp)] [-h|--help] [-ip|--install_path PATH]"
  echo "                    [-j|--num_procs #processes] [-t|test \"TESTOPTS\"] [-v|--verbose]"
  echo ""
  echo "Hints:"
  echo " - Remove the folder ${CP2K_ROOT}/build to (re)build CP2K from scratch"
  echo " - Remove the folder ${CP2K_ROOT}/spack to (re)build CP2K and all its dependencies from scratch (takes a long time)"
  echo " - The folder ${CP2K_ROOT}/install is updated after each successful run"
  echo ""
  ${EXIT_CMD} 0
fi

echo "CMAKE_BUILD_TYPE = ${BUILD_TYPE}"
echo "CP2K_VERSION     = ${CP2K_VERSION}"
echo "INSTALL_PREFIX   = ${INSTALL_PREFIX}"
echo "INSTALL_MESSAGE  = ${INSTALL_MESSAGE}"
echo "NUM_PROCS        = ${NUM_PROCS} (processes)"
echo "RUN_TEST         = ${RUN_TEST}"
if [[ "${RUN_TEST}" == "yes" ]]; then
  echo "TESTOPTS         = \"${TESTOPTS}\""
fi
echo "VERBOSE_FLAG     = ${VERBOSE_FLAG}"
echo "VERBOSE_MAKEFILE = ${VERBOSE_MAKEFILE}"
echo "VERBOSE          = ${VERBOSE}"
if (($# > 0)); then
  echo "Remaining args   =" "$@" "(not used yet)"
fi

# Check if a valid CMake build type is selected
case "${BUILD_TYPE}" in
  Debug | Release | RelWithDebInfo | MinSizeRel)
    true
    ;;
  *)
    echo -e "\nERROR: Invalid CMake build type \"${BUILD_TYPE}\" selected"
    ${EXIT_CMD} 1
    ;;
esac

# Check if all mandatory packages for the Spack/CMake build are installed in the environment
for package in bzip2 g++ gcc gfortran git gmake gzip patch python3 tar wget xz zip; do
  if ! command -v "${package}" &> /dev/null; then
    echo "ERROR: The package \"${package}\" is mandatory to build CP2K with Spack/CMake"
    echo "       Install the missing package and re-run the script"
    ${EXIT_CMD} 1
  fi
done

# Check if all recommended packages for the Spack/CMake build are installed in the environment
for package in autoconf automake bison flex libtool m4 pkg-config; do
  if ! command -v "${package}" &> /dev/null; then
    echo "INFO:  The package \"${package}\" was not found but is recommended with Spack/CMake"
  fi
done

# Check if python3 is available
if ! command -v python3 &> /dev/null; then
  echo "ERROR: python3 NOT found"
  ${EXIT_CMD} 1
fi

# Check if the Pthon version is new enough
if ! python3 -c 'import sys; sys.exit(not(sys.version_info >= (3, 9)))'; then
  echo "ERROR: Python version is NOT >= 3.9 (needed for Spack cache)"
  ${EXIT_CMD} 1
fi

### Build CP2K dependencies with Spack if needed or requested ###

# Spack versions
export SPACK_VERSION="${SPACK_VERSION:-1.1.0}"
export SPACK_PACKAGES_VERSION="${SPACK_PACKAGES_VERSION:-2025.11.0}"

export SPACK_BUILD_PATH="${CP2K_ROOT}/spack"
export SPACK_ROOT="${SPACK_BUILD_PATH}/spack-${SPACK_VERSION}"

# If requested, remove the spack folder for (re)building all CP2K dependencies
[[ "${BUILD_DEPS}" == "yes" ]] && rm -rf "${SPACK_BUILD_PATH}"

if [[ ! -d "${SPACK_BUILD_PATH}" ]]; then

  # Create a new local spack folder
  mkdir -p "${SPACK_BUILD_PATH}"
  cd "${SPACK_BUILD_PATH}" || ${EXIT_CMD} 1

  # The package podman is required for using a local Spack cache
  if command -v podman &> /dev/null; then
    ((VERBOSE > 0)) && echo "Installing virtual environment for Python packages"
    HAS_PODMAN="yes"
    # Create and activate a virtual environment for Python packages
    python3 -m venv "${SPACK_BUILD_PATH}/venv"
    export PATH="${SPACK_BUILD_PATH}/venv/bin:${PATH}"
    pip3 install "${VERBOSE_FLAG}" --upgrade pip
    pip3 install "${VERBOSE_FLAG}" boto3==1.38.11 google-cloud-storage==3.1.0
    # Start Spack cache
    "${CP2K_ROOT}"/tools/docker/spack_cache_start.sh
  else
    HAS_PODMAN="no"
    echo "INFO: podman was not found"
    echo "INFO: Install the package podman to take advantage of a Spack cache"
    echo "INFO: This accelerates a rebuild of the CP2K dependencies, significantly"
  fi

  # Install Spack
  ((VERBOSE > 0)) && echo "Installing Spack ${SPACK_VERSION}"
  export SPACK_REPO="https://github.com/spack/spack"
  if [[ ! -d "${SPACK_ROOT}" ]]; then
    wget -q "${SPACK_REPO}/archive/v${SPACK_VERSION}.tar.gz"
    tar -xzf "v${SPACK_VERSION}.tar.gz"
    rm -f "v${SPACK_VERSION}.tar.gz"
  fi
  export PATH="${SPACK_ROOT}/bin:${PATH}"

  # Isolate user config/cache
  export SPACK_DISABLE_LOCAL_CONFIG=true
  export SPACK_USER_CONFIG_PATH="${SPACK_BUILD_PATH}"
  export SPACK_USER_CACHE_PATH="${SPACK_BUILD_PATH}/cache"
  mkdir -p "${SPACK_USER_CACHE_PATH}"

  # Install Spack packages
  ((VERBOSE > 0)) && echo "Installing Spack packages ${SPACK_PACKAGES_VERSION}"
  export SPACK_PACKAGES_REPO="https://github.com/spack/spack-packages"
  export SPACK_PACKAGES_ROOT="${SPACK_BUILD_PATH}/spack-packages-${SPACK_PACKAGES_VERSION}"
  wget -q "${SPACK_PACKAGES_REPO}/archive/v${SPACK_PACKAGES_VERSION}.tar.gz"
  tar -xzf "v${SPACK_PACKAGES_VERSION}.tar.gz"
  rm -f "v${SPACK_PACKAGES_VERSION}.tar.gz"

  # Add local Spack cache if podman is available
  if [[ "${HAS_PODMAN}" == "yes" ]]; then
    export SPACK_CACHE="s3://spack-cache --s3-endpoint-url=http://localhost:9000"
    if ! spack mirror list | grep -q "local-cache"; then
      echo "Setting up Spack cache"
      "${CP2K_ROOT}"/tools/docker/scripts/setup_spack_cache.sh
    fi
    ((VERBOSE > 0)) && "${CP2K_ROOT}"/tools/docker/spack_cache_list.sh
  else
    echo "INFO: A local Spack cache is NOT used"
  fi

  # Initialize Spack shell hooks
  source "${SPACK_ROOT}/share/spack/setup-env.sh"

  # Prepare the CP2K Spack configuration file
  export CP2K_CONFIG_FILE="${SPACK_BUILD_PATH}/cp2k_deps_${CP2K_VERSION}.yaml"
  sed -e "s|root: /opt/spack|root: ${SPACK_ROOT}/opt/spack|" \
    -e "/\"build_type=/s|build_type=[^\"]*|build_type=${BUILD_TYPE}|" \
    "${CP2K_ROOT}/tools/spack/cp2k_deps_${CP2K_VERSION}.yaml" > "${CP2K_CONFIG_FILE}"

  # Activate CP2K environment if needed
  if spack env list | grep -q "${CP2K_ENV}"; then
    if [[ "${SPACK_ENV:-}" == "${CP2K_ENV}" ]]; then
      echo "The Spack environment ${CP2K_ENV} exists and is activated"
    else
      echo "The Spack environment ${CP2K_ENV} exists but is NOT activated"
      spack env activate "${CP2K_ENV}"
      echo "The Spack environment ${CP2K_ENV} has been activated"
    fi
  else
    echo "The Spack environment ${CP2K_ENV} does NOT exist"
    spack env create "${CP2K_ENV}" "${CP2K_CONFIG_FILE}" && spack env activate "${CP2K_ENV}"
    echo "The Spack environment ${CP2K_ENV} has been created and activated"
  fi

  # Add Spack packages builtin repository when missing
  if ! spack repo list | grep -q "builtin"; then
    spack repo add --scope env:"${CP2K_ENV}" "${SPACK_PACKAGES_ROOT}/repos/spack_repo/builtin"
  else
    echo "Spack built-in repo already present ... skipping add"
  fi

  # Add the local CP2K development Spack repository when missing
  export CP2K_REPO="cp2k_dev"
  if ! spack repo list | grep -q "${CP2K_REPO}"; then
    spack repo add --scope "env:${CP2K_ENV}" "${CP2K_ROOT}/tools/spack/spack_repo/${CP2K_REPO}"
  else
    echo "The CP2K development repo ${CP2K_REPO} is already present ... skipping add"
  fi

  spack -e ${CP2K_ENV} repo list

  # Find all compilers
  spack compiler find

  spack find -c

  # Find all external packages
  spack -C "${SPACK_USER_CONFIG_PATH}" external find --not-buildable

  # Install CP2K dependencies via Spack
  spack -e "${CP2K_ENV}" --no-user-config --no-system-config concretize -f
  spack -e "${CP2K_ENV}" env depfile -o spack_makefile
  make -j"${NUM_PROCS}" --file=spack_makefile SPACK_COLOR=never --output-sync=recurse

  cd "${CP2K_ROOT}" || ${EXIT_CMD} 1

else

  # Initialize Spack shell hooks
  source "${SPACK_ROOT}"/share/spack/setup-env.sh

  if [[ "${HAS_PODMAN}" == "yes" ]]; then
    # Start Spack cache
    "${CP2K_ROOT}"/tools/docker/spack_cache_start.sh
    podman container list
    ((VERBOSE > 0)) && "${CP2K_ROOT}"/tools/docker/spack_cache_list.sh
  fi

fi

### End of build CP2K dependencies ###

### Build CP2K ###

spack env list
spack env status
eval "$(spack env activate --sh ${CP2K_ENV})"

# CMake configuration step
export CMAKE_BUILD_PATH="${CP2K_ROOT}/build"

# PyTorch's TorchConfig.cmake is buried in the Python site-packages directory
Torch_DIR="$(dirname "$(find "${SPACK_ROOT}" ! -type l -name TorchConfig.cmake | tail -n 1)")"
export Torch_DIR

if [[ ! -d "${CMAKE_BUILD_PATH}" ]]; then
  mkdir -p "${CMAKE_BUILD_PATH}"
  case "${CP2K_VERSION}" in
    "psmp")
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" \
        -GNinja \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH=ON \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        -DCP2K_USE_EVERYTHING=ON \
        -DCP2K_USE_DLAF=OFF \
        -DCP2K_USE_LIBXSMM=ON \
        -DCP2K_USE_MPI=ON \
        -DCP2K_USE_OPENPMD=ON \
        -DCP2K_USE_TBLITE=OFF \
        -Werror=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    "ssmp")
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" \
        -GNinja \
        -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH=ON \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        -DCP2K_USE_EVERYTHING=ON \
        -DCP2K_USE_DLAF=OFF \
        -DCP2K_USE_LIBXSMM=ON \
        -DCP2K_USE_MPI=OFF \
        -DCP2K_USE_OPENPMD=OFF \
        -DCP2K_USE_TBLITE=OFF \
        -Werror=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    *)
      echo "ERROR: Unknown CP2K version ${CP2K_VERSION} specified"
      ${EXIT_CMD} 1
      ;;
  esac
fi

if ((EXIT_CODE != 0)); then
  echo "ERROR: The CMake configuration step failed with the error code ${EXIT_CODE}"
  echo "       You can try to remove the build folder with 'rm -rf build' and re-run"
  ${EXIT_CMD} "${EXIT_CODE}"
fi

# CMake build step
echo -e '\n*** Compiling CP2K ***\n'
cmake --build "${CMAKE_BUILD_PATH}" --parallel "${NUM_PROCS}" -- "${VERBOSE_FLAG}" |& tee "${CMAKE_BUILD_PATH}/ninja.log"
EXIT_CODE=$?
if ((EXIT_CODE != 0)); then
  echo "ERROR: The CMake build step failed with the error code ${EXIT_CODE}"
  ${EXIT_CMD} "${EXIT_CODE}"
fi

# CMake install step
echo -e '\n*** Installing CP2K ***\n'
cmake --install "${CMAKE_BUILD_PATH}" |& tee "${CMAKE_BUILD_PATH}/install.log"
EXIT_CODE=$?
if ((EXIT_CODE != 0)); then
  echo "ERROR: The CMake installation step failed with the error code ${EXIT_CODE}"
  ${EXIT_CMD} "${EXIT_CODE}"
fi

# Retrieve paths to "hidden" libraries
if ldd "${INSTALL_PREFIX}/bin/cp2k.${CP2K_VERSION}" | grep -q "not found"; then
  echo -e "\n*** Some libraries referenced by the CP2K binary are NOT found:"
  ldd "${INSTALL_PREFIX}/bin/cp2k.${CP2K_VERSION}" | grep "not found"
  echo -e "\n*** Trying to update the LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}\n"
  for libname in $(ldd "${INSTALL_PREFIX}/bin/cp2k.${CP2K_VERSION}" | grep "not found" | awk '{print $1}' | sort | uniq); do
    library=$(find -L "${SPACK_ROOT}/opt/spack/view" "${INSTALL_PREFIX}/lib" -name "${libname}" | tail -n 1)
    library_path="$(dirname "${library}")"
    if [[ -n "${library_path}" ]]; then
      case ":${LD_LIBRARY_PATH}:" in
        *":${library_path}:"*)
          echo "The library path ${library_path} is already present ... skipping"
          ;;
        *)
          LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}${library_path}"
          echo "The library path ${library_path} was appended to LD_LIBRARY_PATH"
          ;;
      esac
    fi
  done
  if ! ldd "${INSTALL_PREFIX}/bin/cp2k.${CP2K_VERSION}" | grep -q "not found"; then
    echo -e "\n*** The update of the LD_LIBRARY_PATH was successful"
    echo -e "\n*** New LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"
  else
    echo -e "\n*** The update of the LD_LIBRARY_PATH was NOT successful"
    echo -e "\n*** The following libraries were NOT found:"
    ldd "${INSTALL_PREFIX}/bin/cp2k.${CP2K_VERSION}" | grep "not found" | sort | uniq
  fi
  export LD_LIBRARY_PATH
fi

# Prepare script to run the regression tests
TESTOPTS="--cp2kdatadir ${INSTALL_PREFIX}/share/cp2k/data  --maxtasks ${NUM_PROCS} --workbasedir ${INSTALL_PREFIX}/regtesting ${TESTOPTS}"
cat << *** > "${INSTALL_PREFIX}"/bin/run_tests
#!/bin/bash
ulimit -c 0 -s unlimited
export OMP_STACKSIZE=64M
export PATH=${INSTALL_PREFIX}/bin:${PATH}
export LD_LIBRARY_PATH=${INSTALL_PREFIX}/lib:${LD_LIBRARY_PATH}
ldd "${INSTALL_PREFIX}/bin/cp2k.${CP2K_VERSION}" | grep "not found" | sort | uniq
${CP2K_ROOT}/tests/do_regtest.py ${TESTOPTS} \$* ${INSTALL_PREFIX}/bin ${CP2K_VERSION}
***
chmod 750 "${INSTALL_PREFIX}/bin/run_tests"

# Optionally, launch test run
if [[ "${RUN_TEST}" == "yes" ]]; then
  echo -e "\n*** Launching regression test run using the script ${INSTALL_PREFIX}/bin/run_tests\n"
  "${INSTALL_PREFIX}"/bin/run_tests
  EXIT_CODE=$?
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The regression test run failed with the error code ${EXIT_CODE}"
    ${EXIT_CMD} "${EXIT_CODE}"
  fi
else
  echo -e "\n*** A regression test run can be launched using the script ${INSTALL_PREFIX}/bin/run_tests\n"
fi
