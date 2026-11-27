#!/usr/bin/env bash
#
# build_cp2k.sh: CMake configuration, build, and installation of CP2K
#
# Usage: build_cp2k.sh CONFIG_FILE...
#        build_cp2k.sh [OPTIONS]
#
# This script is used for both the Spack-based and the toolchain-based CP2K builds:
#
#  - With configuration files (options.conf and deps.conf written by make_cp2k.sh
#    and build_deps.sh), the settings are taken from these files and the Spack
#    environment with the CP2K dependencies is activated (step 3 of make_cp2k.sh).
#    After the installation, the environment of the CP2K binaries is set up
#    (launch and run_tests scripts) and the settings provided by the CP2K build
#    are written to CONFIG_PATH/build.conf.
#  - Without configuration files, CP2K is only configured, built, and installed
#    with the settings given by the flags, e.g. for a toolchain build:
#
#      build_cp2k.sh --setup /opt/cp2k-toolchain/install/setup \
#                    --preset native-gnu-x86_64 --install_prefix /opt/cp2k \
#                    --cmake_args "-DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF"
#
# The flags cannot be combined with configuration files.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

CP2K_ROOT_DIR="$(realpath "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/../..")"

# shellcheck source=/dev/null
source "${CP2K_ROOT_DIR}/tools/make_cp2k/common.sh"

print_usage() {
  echo "Usage: $(basename "${0}") CONFIG_FILE..."
  echo "       $(basename "${0}") [OPTIONS]"
  echo ""
  echo "Options (only without configuration files):"
  echo " -s, --setup FILE          : Source FILE to set up the environment, e.g. the toolchain setup file"
  echo " -p, --preset NAME         : CMake configure preset (default: none)"
  echo " -c, --cmake_args ARGS     : Further CMake arguments as one string, e.g. \"-DCP2K_USE_MPI=OFF -DCP2K_USE_LIBXC=ON\""
  echo "                             (\${VAR} references are expanded after sourcing the setup file)"
  echo " -b, --build_dir DIR       : CMake build folder (default: ${CP2K_ROOT_DIR}/build)"
  echo " -i, --install_prefix DIR  : CMake installation prefix (default: ${CP2K_ROOT_DIR}/install)"
  echo " -j, --jobs N              : Number of parallel build jobs (default: $(nproc))"
  echo " -fw, --fail_on_warnings   : Fail if the CMake configuration step issues any warning"
  echo " -ni, --no_install         : Skip the installation step and only link the data folder to"
  echo "                             <install_prefix>/share/cp2k/data (binaries remain in <build_dir>/bin)"
  echo " -q, --quiet               : Write the build output only to <build_dir>/ninja.log, print the end"
  echo "                             of the log and copy it to /workspace/artifacts if the build fails,"
  echo "                             and report failures with Summary: and Status: lines (for the CP2K CI)"
  echo " -v, --verbose             : Print the full compiler and linker command lines"
  echo " -h, --help                : Print this help information"
}

if (($# == 0)); then
  echo "ERROR: Neither configuration files nor flags specified"
  print_usage
  exit 1
fi

# Parse the flags
SETUP_FILE=""
CMAKE_PRESET_FLAG=""
CMAKE_ARGS=""
BUILD_DIR=""
INSTALL_DIR=""
JOBS=""
FAIL_ON_WARNINGS="no"
NO_INSTALL="no"
QUIET="no"
VERBOSE_BUILD="no"
FLAGS_FOUND="no"
CONFIG_FILES=()
while (($# > 0)); do
  case "${1}" in
    -s | --setup | -p | --preset | -c | --cmake_args | -b | --build_dir | -i | --install_prefix | \
      -j | --jobs)
      if (($# < 2)) || [[ -z "${2}" ]]; then
        echo "ERROR: No argument found for flag \"${1}\""
        exit 1
      fi
      case "${1}" in
        -s | --setup) SETUP_FILE="${2}" ;;
        -p | --preset) CMAKE_PRESET_FLAG="${2}" ;;
        -c | --cmake_args) CMAKE_ARGS="${2}" ;;
        -b | --build_dir) BUILD_DIR="${2}" ;;
        -i | --install_prefix) INSTALL_DIR="${2}" ;;
        -j | --jobs) JOBS="${2}" ;;
      esac
      FLAGS_FOUND="yes"
      shift 2
      ;;
    -fw | --fail_on_warnings | -ni | --no_install | -q | --quiet | -v | --verbose)
      case "${1}" in
        -fw | --fail_on_warnings) FAIL_ON_WARNINGS="yes" ;;
        -ni | --no_install) NO_INSTALL="yes" ;;
        -q | --quiet) QUIET="yes" ;;
        -v | --verbose) VERBOSE_BUILD="yes" ;;
      esac
      FLAGS_FOUND="yes"
      shift 1
      ;;
    -h | --help)
      print_usage
      exit 0
      ;;
    -*)
      echo "ERROR: Unknown flag \"${1}\" found"
      print_usage
      exit 1
      ;;
    *)
      CONFIG_FILES+=("${1}")
      shift 1
      ;;
  esac
done

# Print an error message and also the Summary and Status lines for the CP2K CI if requested
report_failure() {
  echo -e "\nERROR: ${1}"
  if [[ "${QUIET}" == "yes" ]]; then
    echo -e "\nSummary: ${2}"
    echo -e "Status: FAILED\n"
  fi
}

# Expand the ${VAR} references in the given string without evaluating anything else
expand_variables() {
  local string="${1}" name
  while [[ "${string}" =~ \$\{([A-Za-z_][A-Za-z0-9_]*)\} ]]; do
    name="${BASH_REMATCH[1]}"
    if [[ ! -v ${name} ]]; then
      echo "ERROR: The variable ${name} referenced in \"${1}\" is not set" >&2
      return 1
    fi
    string="${string//"\${${name}}"/${!name}}"
  done
  echo "${string}"
}

if ((${#CONFIG_FILES[@]} > 0)); then
  if [[ "${FLAGS_FOUND}" == "yes" ]]; then
    echo "ERROR: The flags cannot be combined with configuration files"
    print_usage
    exit 1
  fi
  load_config "${CONFIG_FILES[@]}"
  activate_spack_env
  spack env status
  set_cmake_werror_flags
  NINJA_OUTPUT_FLAG="${VERBOSE_FLAG}"
else
  # Set up the environment for a build without configuration files
  if [[ -n "${SETUP_FILE}" ]]; then
    if [[ ! -f "${SETUP_FILE}" ]]; then
      echo "ERROR: The setup file ${SETUP_FILE} was not found"
      exit 1
    fi
    # The setup file is not necessarily safe for unset variables
    set +u
    # shellcheck source=/dev/null
    source "${SETUP_FILE}"
    set -u
  fi
  if [[ -n "${JOBS}" ]] && [[ ! "${JOBS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: Invalid number of build jobs \"${JOBS}\" (flag --jobs)"
    exit 1
  fi
  CMAKE_ARGS="$(expand_variables "${CMAKE_ARGS}")" || exit 1
  read -ra CMAKE_ARGV <<< "${CMAKE_ARGS}"
  if [[ -n "${CMAKE_PRESET_FLAG}" ]]; then
    CMAKE_ARGV=(--preset "${CMAKE_PRESET_FLAG}" "${CMAKE_ARGV[@]}")
  fi
  # The variable names used by the build and installation steps below
  CP2K_ROOT="${CP2K_ROOT_DIR}"
  CMAKE_BUILD_PATH="${BUILD_DIR:-${CP2K_ROOT_DIR}/build}"
  INSTALL_PREFIX="${INSTALL_DIR:-${CP2K_ROOT_DIR}/install}"
  NUM_PROCS="${JOBS:-$(nproc)}"
  if [[ "${VERBOSE_BUILD}" == "yes" ]]; then
    NINJA_OUTPUT_FLAG="--verbose"
  else
    NINJA_OUTPUT_FLAG="--quiet"
  fi
  # Record the CP2K revision when building from a source tree without git metadata
  if [[ -n "${GIT_COMMIT_SHA:-}" ]]; then
    echo "git:${GIT_COMMIT_SHA::7}" > "${CP2K_ROOT_DIR}"/REVISION
  fi
fi

### Build CP2K ###

if ((${#CONFIG_FILES[@]} == 0)); then
  # CMake configuration step from scratch for a build without configuration files
  echo -e "\n*** Configuring CP2K ***\n"
  mkdir -p "${CMAKE_BUILD_PATH}"
  cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" --fresh \
    -GNinja \
    -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
    "${CMAKE_ARGV[@]}" |&
    tee "${CMAKE_BUILD_PATH}/cmake.log"
  EXIT_CODE=${PIPESTATUS[0]}
  if ((EXIT_CODE != 0)); then
    report_failure "The CMake configuration step failed with the error code ${EXIT_CODE}" "CMake failed"
    exit "${EXIT_CODE}"
  fi
  if [[ "${FAIL_ON_WARNINGS}" == "yes" ]] && grep -A5 'CMake Warning' "${CMAKE_BUILD_PATH}/cmake.log"; then
    report_failure "The CMake configuration step issued warnings" "Found CMake warnings"
    exit 1
  fi
elif [[ ! -d "${CMAKE_BUILD_PATH}" ]]; then
  mkdir -p "${CMAKE_BUILD_PATH}"
  case "${CP2K_VERSION}" in
    pdbg | psmp | sdbg | ssmp)
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" --preset "${CMAKE_PRESET}" \
        -GNinja \
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} \
        -DCMAKE_BUILD_TYPE="${CP2K_BUILD_TYPE}" \
        -DCMAKE_Fortran_COMPILER_LAUNCHER="${Fortran_COMPILER_LAUNCHER}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH="ON" \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        ${CMAKE_FEATURE_FLAGS} \
        ${CMAKE_CUDA_FLAGS} \
        "${CMAKE_WERROR_FLAGS[@]}" |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    ssmp-static)
      # Find some static libraries in advance
      LIBOPENBLAS=$(find -L "${SPACK_ROOT}"/opt/spack/view -name libopenblas.a)
      OPENBLAS_INCLUDE_DIR="$(
        dirname "$(find -L "${SPACK_ROOT}"/opt/spack/view -name cblas.h -print -quit)"
      )"
      LIBM="$(find /usr -name libm.a 2> /dev/null)"
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" --preset "${CMAKE_PRESET}" \
        -GNinja \
        -DBUILD_SHARED_LIBS="OFF" \
        -DCMAKE_BUILD_TYPE="${CP2K_BUILD_TYPE}" \
        -DCMAKE_EXE_LINKER_FLAGS="-static" \
        -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH="ON" \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        -DCP2K_BLAS_VENDOR="OpenBLAS" \
        -DCP2K_BLAS_INCLUDE_DIRS="${OPENBLAS_INCLUDE_DIR}" \
        -DCP2K_BLAS_LINK_LIBRARIES="${LIBOPENBLAS};${LIBM}" \
        -DCP2K_LAPACK_LINK_LIBRARIES="${LIBOPENBLAS};${LIBM}" \
        ${CMAKE_FEATURE_FLAGS} \
        "${CMAKE_WERROR_FLAGS[@]}" |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      # It is almost impossible to avoid that shared libraries are pulled in
      # when the CP2K dependencies are built and thus just change all shared
      # to static libraries hoping that these are also available
      if ((EXIT_CODE == 0)); then
        sed -E -e "s/\.so/.a/g" -i "${CMAKE_BUILD_PATH}"/build.ninja
      fi
      ;;
    *)
      echo "ERROR: Unknown CP2K version ${CP2K_VERSION} specified"
      exit 1
      ;;
  esac
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The CMake configuration step failed with the error code ${EXIT_CODE}"
    echo "       You can try to remove the build folder with 'rm -rf build' and re-run"
    echo "       or even start the whole CP2K installation from scratch with the \"-bd\" flag"
    exit "${EXIT_CODE}"
  fi
fi

# CMake build step
echo -e "\n*** Compiling CP2K ***\n"
if [[ "${QUIET}" == "yes" ]]; then
  cmake --build "${CMAKE_BUILD_PATH}" --parallel "${NUM_PROCS}" -- "${NINJA_OUTPUT_FLAG}" &> "${CMAKE_BUILD_PATH}"/ninja.log
else
  cmake --build "${CMAKE_BUILD_PATH}" --parallel "${NUM_PROCS}" -- "${NINJA_OUTPUT_FLAG}" |& tee "${CMAKE_BUILD_PATH}"/ninja.log
fi
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  if [[ "${QUIET}" == "yes" ]]; then
    # Stay below the 3 MB limit for the CP2K dashboard reports
    tail -n 500 "${CMAKE_BUILD_PATH}"/ninja.log
    mkdir -p /workspace/artifacts
    cp "${CMAKE_BUILD_PATH}"/ninja.log /workspace/artifacts/
  fi
  report_failure "The CMake build step failed with the error code ${EXIT_CODE}" "Compilation failed"
  exit "${EXIT_CODE}"
fi

# Only link the data folder to the installation prefix if the installation is skipped,
# since CP2K looks for its data files in <install_prefix>/share/cp2k/data by default
if [[ "${NO_INSTALL}" == "yes" ]]; then
  if [[ ! -e "${INSTALL_PREFIX}"/share/cp2k/data ]]; then
    mkdir -p "${INSTALL_PREFIX}"/share/cp2k
    ln -s "${CP2K_ROOT}"/data "${INSTALL_PREFIX}"/share/cp2k/data
  fi
  echo -e "\n*** CP2K was built in ${CMAKE_BUILD_PATH} without installation ***\n"
  exit 0
fi

# CMake install step
echo -e "\n*** Installing CP2K ***\n"
cmake --install "${CMAKE_BUILD_PATH}" |& tee "${CMAKE_BUILD_PATH}"/install.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  report_failure "The CMake installation step failed with the error code ${EXIT_CODE}" "Installation failed"
  exit "${EXIT_CODE}"
fi

# The following set-up of the environment for the CP2K binaries requires the configuration files
if ((${#CONFIG_FILES[@]} == 0)); then
  echo -e "\n*** CP2K was installed in ${INSTALL_PREFIX} ***\n"
  exit 0
fi

# Install Skala resource files if needed
export SKALA_MODEL="(not available)"
unset SKALA_CUDA_MODEL
if spack location -i gauxc &> /dev/null; then
  GAUXC_PATH="share/gauxc/onedft_models"
  GAUXC_PREFIX="$(spack location -i gauxc)"
  GAUXC_MODEL_SOURCE_PATH="${GAUXC_PREFIX%/}/${GAUXC_PATH}"
  if [[ -d "${GAUXC_MODEL_SOURCE_PATH}" ]]; then
    GAUXC_MODEL_TARGET_PATH="${INSTALL_PREFIX%/}/${GAUXC_PATH}"
    mkdir -p "${GAUXC_MODEL_TARGET_PATH}"
    if compgen -G "${GAUXC_MODEL_SOURCE_PATH}/*.fun" > /dev/null; then
      cp "${GAUXC_MODEL_SOURCE_PATH}"/*.fun "${GAUXC_MODEL_TARGET_PATH}/"
    else
      echo -e "\nERROR: No GauXC model files found in source folder ${GAUXC_MODEL_SOURCE_PATH}"
      exit 1
    fi
    # Distinguish the CPU model (skala*.fun) from the CUDA model (skala*-cuda.fun).
    # Both match a plain skala* glob, so assign them to the dedicated variables
    # that CP2K expects: SKALA_MODEL for host and SKALA_CUDA_MODEL for device runs.
    shopt -s nullglob
    SKALA_MODEL_FILES=("${GAUXC_MODEL_TARGET_PATH}"/skala*.fun)
    shopt -u nullglob
    for skala_file in "${SKALA_MODEL_FILES[@]}"; do
      case "$(basename "${skala_file}")" in
        *-cuda.fun)
          export SKALA_CUDA_MODEL="${skala_file}"
          ;;
        *)
          export SKALA_MODEL="${skala_file}"
          ;;
      esac
    done
    if [[ "${SKALA_MODEL}" == "(not available)" ]]; then
      echo -e "\nERROR: Failed to resolve SKALA_MODEL in target path ${GAUXC_MODEL_TARGET_PATH}"
      exit 1
    fi
    if ((CUDA_SM_CODE > 0)) && [[ -z "${SKALA_CUDA_MODEL:-}" ]]; then
      echo -e "\nERROR: CUDA support is enabled but no SKALA_CUDA_MODEL file was found in ${GAUXC_MODEL_TARGET_PATH}"
      exit 1
    fi
  fi
fi
echo -e "\nSKALA_MODEL = ${SKALA_MODEL}"
echo -e "SKALA_CUDA_MODEL = ${SKALA_CUDA_MODEL:-not available}"

# Emit the Skala model exports for the generated environment scripts
write_skala_model_exports() {
  echo "export SKALA_MODEL=\"${SKALA_MODEL}\""
  if [[ -n "${SKALA_CUDA_MODEL:-}" ]]; then
    echo "export SKALA_CUDA_MODEL=\"${SKALA_CUDA_MODEL}\""
  fi
}

# Collect and compress all log files when building within a container
if [[ "${IN_CONTAINER}" == "yes" ]]; then
  if ! cat "${CMAKE_BUILD_PATH}"/cmake.log \
    "${CMAKE_BUILD_PATH}"/ninja.log \
    "${CMAKE_BUILD_PATH}"/install.log |
    gzip > "${CP2K_ROOT}"/build_cp2k.log.gz; then
    echo -e "\nERROR: The compressed log file generation failed"
    exit 1
  fi
fi

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
cd "${INSTALL_PREFIX}"/bin || exit 1
for binary in *."${VERSION}"; do
  if ! ln -sf "${binary}" "${binary%%."${VERSION}"}"; then
    echo "ERROR: The creation of a symbolic link for the binary ${binary} failed"
    exit
  fi
done
ln -sf cp2k."${VERSION}" cp2k_shell
cd "${CP2K_ROOT}" || exit 1

# Configure Open MPI environment for local CP2K runs
if [[ "${MPI_MODE}" == "openmpi" ]]; then
  OMPI_VARS="export OMPI_MCA_plm_rsh_agent=/bin/false PRTE_MCA_hwloc_default_binding_policy=none"
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    OMPI_VARS="${OMPI_VARS} OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1"
    OMPI_VARS="${OMPI_VARS} OMPI_MCA_mpi_yield_when_idle=1 OMPI_MCA_btl=self,sm OMPI_MCA_pml=ob1"
  fi
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

# Install LSAN suppressions (only needed for dbg and ASAN binaries)
cp "${CP2K_ROOT}"/tools/spack/lsan.supp "${INSTALL_PREFIX}"/bin

# Assemble flags for running the regression tests
TESTOPTS="--cp2kdatadir ${INSTALL_PREFIX}/share/cp2k/data  --maxtasks ${NUM_PROCS} --workbasedir ${INSTALL_PREFIX}/regtesting ${TESTOPTS}"
export TESTOPTS

# Create launch script for environment setup
LAUNCH_SCRIPT="${INSTALL_PREFIX}"/bin/launch
export LAUNCH_SCRIPT
{
  cat << ***
#!/bin/bash
ulimit -c 0 -s unlimited
export ASAN_OPTIONS="detect_leaks=1"
export LSAN_OPTIONS="suppressions=${INSTALL_PREFIX}/bin/lsan.supp"
export PATH="${INSTALL_PREFIX}/bin:${INSTALL_PREFIX}/ase/bin:${PATH}"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-2}
export OMP_STACKSIZE=256M
[[ -f ${INSTALL_PREFIX}/ase/config.ini ]] && export ASE_CONFIG_PATH="${INSTALL_PREFIX}/ase/config.ini"
[[ -f ${INSTALL_PREFIX}/bin/GMXRC ]] && source ${INSTALL_PREFIX}/bin/GMXRC
${OMPI_VARS}
***
  write_skala_model_exports
  cat << ***
exec "\$@"
***
} > "${LAUNCH_SCRIPT}"
chmod 750 "${LAUNCH_SCRIPT}"

# Create shortcut for launching the CP2K regression tests
{
  cat << ***
#!/bin/bash
if [[ "${VERSION}" =~ ^(s|p)dbg$ ]] || [[ "${CMAKE_PRESET}" == "asan" ]]; then
  echo "ASAN_OPTIONS = \${ASAN_OPTIONS}"
  echo "LSAN_OPTIONS = \${LSAN_OPTIONS}"
fi
ldd -- ${INSTALL_PREFIX}/bin/cp2k.${VERSION} 2>&1 | grep -E 'not ' | sort | uniq
***
  write_skala_model_exports
  cat << ***
${CP2K_ROOT}/tests/do_regtest.py ${TESTOPTS} \$* ${INSTALL_PREFIX}/bin ${VERSION}
***
} > "${INSTALL_PREFIX}"/bin/run_tests
chmod 750 "${INSTALL_PREFIX}"/bin/run_tests

# Write the settings provided by the CP2K build
write_config LAUNCH_SCRIPT LD_LIBRARY_PATH SKALA_MODEL SKALA_CUDA_MODEL > "${CONFIG_PATH}/build.conf"

# Print usage examples if no regression test run is requested
if [[ "${RUN_TEST}" != "yes" ]]; then
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    if ((CUDA_SM_CODE > 0)); then
      DEVICE_FLAG=" --device nvidia.com/gpu=all"
    else
      DEVICE_FLAG=""
    fi
    echo ""
    echo "*** A regression test run can be launched with"
    echo "    podman run -it${DEVICE_FLAG} --rm ${IMAGE_TAG} run_tests"
    echo ""
    if [[ "${VERSION}" == @(sdbg|ssmp) ]]; then
      echo "*** A CP2K run using 2 OpenMP threads (default) can be launched with"
      echo "    podman run -it${DEVICE_FLAG} --rm ${IMAGE_TAG} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    else
      echo "*** An MPI-parallel CP2K run using 2 OpenMP threads for each of the 4 MPI ranks can be launched with"
      echo "    podman run -it${DEVICE_FLAG} --rm ${IMAGE_TAG} mpiexec -n 4 ${ENV_VAR_FLAG} OMP_NUM_THREADS=2 cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    fi
  else
    echo ""
    echo "*** A regression test run can be launched with"
    echo "    ${LAUNCH_SCRIPT} run_tests"
    echo ""
    if [[ "${VERSION}" == @(sdbg|ssmp) ]]; then
      echo "*** A CP2K run using 2 OpenMP threads (default) can be launched with"
      echo "    ${LAUNCH_SCRIPT} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
      echo ""
      echo "*** A CP2K run using only 4 OpenMP threads can be launched with"
      echo "    export OMP_NUM_THREADS=4; ${LAUNCH_SCRIPT} cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    else
      echo "*** An MPI-parallel CP2K run using 2 OpenMP threads for each of the 4 MPI ranks can be launched with"
      echo "    export OMP_NUM_THREADS=2; ${LAUNCH_SCRIPT} mpiexec -n 4 cp2k ${CP2K_ROOT}/benchmarks/CI/H2O-32_md.inp"
    fi
  fi
fi
