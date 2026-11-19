#!/usr/bin/env bash
#
# build_cp2k.sh: CMake configuration, build, and installation of CP2K (step 3 of make_cp2k.sh)
#
# Usage: build_cp2k.sh CONFIG_FILE...
#
# The settings are read from the configuration files options.conf and deps.conf
# written by make_cp2k.sh and build_deps.sh. The settings provided by the CP2K
# build are written to CONFIG_PATH/build.conf.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"
activate_spack_env
spack env status

### Build CP2K ###

if [[ ! -d "${CMAKE_BUILD_PATH}" ]]; then
  mkdir -p "${CMAKE_BUILD_PATH}"
  case "${CP2K_VERSION}" in
    pdbg | psmp)
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
        -Wno-error=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    sdbg | ssmp)
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
        -Wno-error=dev |&
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
        -Wno-error=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      # It is almost impossible to avoid that shared libraries are pulled in
      # when the CP2K dependencies are built and thus just change all shared
      # to static libraries hoping that these are also available
      sed -E -e "s/\.so/.a/g" -i build/build.ninja
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
cmake --build "${CMAKE_BUILD_PATH}" --parallel "${NUM_PROCS}" -- "${VERBOSE_FLAG}" |& tee "${CMAKE_BUILD_PATH}"/ninja.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  echo "ERROR: The CMake build step failed with the error code ${EXIT_CODE}"
  exit "${EXIT_CODE}"
fi

# CMake install step
echo -e "\n*** Installing CP2K ***\n"
cmake --install "${CMAKE_BUILD_PATH}" |& tee "${CMAKE_BUILD_PATH}"/install.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  echo -e "\nERROR: The CMake installation step failed with the error code ${EXIT_CODE}"
  exit "${EXIT_CODE}"
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
