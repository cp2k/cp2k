# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Runtime setup: library paths, launch and run_tests scripts, benchmark script
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

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
cat << *** > "${LAUNCH_SCRIPT}"
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
export SKALA_MODEL=${SKALA_MODEL}
exec "\$@"
***
chmod 750 "${LAUNCH_SCRIPT}"

# Create shortcut for launching the CP2K regression tests
cat << *** > "${INSTALL_PREFIX}"/bin/run_tests
#!/bin/bash
if [[ "${VERSION}" =~ ^(s|p)dbg$ ]] || [[ "${CMAKE_PRESET}" == "asan" ]]; then
  echo "ASAN_OPTIONS = \${ASAN_OPTIONS}"
  echo "LSAN_OPTIONS = \${LSAN_OPTIONS}"
fi
ldd -- ${INSTALL_PREFIX}/bin/cp2k.${VERSION} 2>&1 | grep -E 'not ' | sort | uniq
export SKALA_MODEL=${SKALA_MODEL}
${CP2K_ROOT}/tests/do_regtest.py ${TESTOPTS} \$* ${INSTALL_PREFIX}/bin ${VERSION}
***
chmod 750 "${INSTALL_PREFIX}"/bin/run_tests

# Collect information from coding convention checks
if [[ "${CHECK_CONVENTIONS}" == "yes" ]]; then
  "${CP2K_ROOT}"/tools/conventions/analyze_gfortran_ast.py "${CMAKE_BUILD_PATH}"/*.ast --jobs "${NUM_PROCS}" \
    &> "${CMAKE_BUILD_PATH}"/ast.issues
  ((VERBOSE > 0)) && cat "${CMAKE_BUILD_PATH}"/ast.issues
  "${CP2K_ROOT}"/tools/conventions/analyze_gfortran_warnings.py "${CMAKE_BUILD_PATH}"/*.warn --jobs "${NUM_PROCS}" \
    &> "${CMAKE_BUILD_PATH}"/warn.issues
  ((VERBOSE > 0)) && cat "${CMAKE_BUILD_PATH}"/warn.issues
  "${CP2K_ROOT}"/tools/conventions/summarize_issues.py --suppressions="${CP2K_ROOT}/tools/conventions/conventions.supp" "${CMAKE_BUILD_PATH}"/*.issues
  cat << *** > "${INSTALL_PREFIX}"/bin/summarize_issues
#!/bin/bash
${CP2K_ROOT}/tools/conventions/summarize_issues.py --suppressions=${CP2K_ROOT}/tools/conventions/conventions.supp ${CMAKE_BUILD_PATH}/*.issues
***
  chmod 750 "${INSTALL_PREFIX}"/bin/summarize_issues
fi

# Create script to run the CP2K benchmarks for psmp builds
if [[ "${VERSION}" == "psmp" ]]; then
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    BENCHMARK_OUTPUT_DIR="/workspace/artifacts"
  else
    BENCHMARK_OUTPUT_DIR="${INSTALL_PREFIX}"/performance_tests
  fi
  cat << *** > "${INSTALL_PREFIX}"/bin/run_benchmarks
#!/bin/bash -e

BENCHMARK_PROFILE=\${1:-openmp}

MAX_PROCS=\$(nproc --all)
if ((MAX_PROCS < 32)); then
  echo -e "\nWARNING: 32 CPU cores are needed for a CP2K benchmark run but if seems only \${MAX_PROCS} are available"
fi

run_benchmark_input() {
  set +e
  local omp_threads="\$1"
  local mpi_ranks="\$2"
  local input="\$3"
  local output="\$4"
  printf "Running %s with %s threads and %s ranks... " \
    "\${input}" "\${omp_threads}" "\${mpi_ranks}"
  if OMP_NUM_THREADS="\${omp_threads}" mpiexec -n "\${mpi_ranks}" cp2k.psmp "\${input}" >"\${output}" 2>&1; then
    echo "done."
  else
    echo "failed."
    echo
    [[ -f "\${output}" ]] && tail -n 1000 "\${output}"
    echo
    echo "Summary: Running \${input} failed."
    echo "Status: FAILED"
    exit 0
  fi
  set -e
}

echo -e "\n========== Running Performance Test =========="
OUTPUT_DIR=${BENCHMARK_OUTPUT_DIR}
mkdir -p \${OUTPUT_DIR}
TIME_START=\$(date +%s)

BENCHMARKS=(
  "${CP2K_ROOT}/benchmarks/QS/H2O-64.inp"
  "${CP2K_ROOT}/benchmarks/QS/H2O-64_nonortho.inp"
  "${CP2K_ROOT}/benchmarks/QS_reference/w64PBE.inp"
  "${CP2K_ROOT}/benchmarks/QS_reference/w64GAPW-SVP.inp"
  "${CP2K_ROOT}/benchmarks/QS_reference/w64SCAN.inp"
  "${CP2K_ROOT}/benchmarks/QS_kp/ZnO.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/H2O-hyb.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/GW_PBE_4benzene.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/RI-HFX_H2O-32.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/RI-MP2_ammonia.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/diag_cu144_broy.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/bench_dftb.inp"
  "${CP2K_ROOT}/benchmarks/QS_single_node/dbcsr.inp"
  "${CP2K_ROOT}/benchmarks/QMMM/MQAE/MQAE_single_node.inp"
)

if [[ "\${BENCHMARK_PROFILE}" == "openmp" ]]; then
  echo 'Plot: name="total_timings_32omp", title="Total Timings with 32 OpenMP Threads", ylabel="time [s]"'
  echo 'Plot: name="total_timings_32mpi", title="Total Timings with 32 MPI Ranks", ylabel="time [s]"'
  echo ""

  for INPUT in "\${BENCHMARKS[@]}"; do
    INPUT_BASENAME=\$(basename "\${INPUT}")
    LABEL=\${INPUT_BASENAME%.*}
    OUTPUT_MPI="\${OUTPUT_DIR}/\${LABEL}_32mpi.out"
    OUTPUT_OMP="\${OUTPUT_DIR}/\${LABEL}_32omp.out"
    cd "\$(dirname "\${INPUT}")"
    run_benchmark_input 1 32 "\${INPUT_BASENAME}" "\${OUTPUT_MPI}"
    run_benchmark_input 32 1 "\${INPUT_BASENAME}" "\${OUTPUT_OMP}"
    cd ..
    echo ""
    ${CP2K_ROOT}/tools/docker/scripts/plot_performance.py \
      "\${LABEL} with 32 OpenMP Threads" "\${LABEL}" "32omp" "\${OUTPUT_OMP}" \
      "\${LABEL} with 32 MPI Ranks" "\${LABEL}" "32mpi" "\${OUTPUT_MPI}"
    echo ""
  done

else
  echo "ERROR: Unknown benchmark profile \${BENCHMARK_PROFILE} specified"
  exit 1
fi

TIME_END=\$(date +%s)
DURATION=\$(printf "%i" \$(((TIME_END - TIME_START) / 60)))

echo -e "\nSummary: Performance test took \${DURATION} minutes."
echo -e "Status: OK\n"
***
  chmod 750 "${INSTALL_PREFIX}"/bin/run_benchmarks
fi
