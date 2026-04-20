#!/bin/bash -e

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
TOOLCHAIN_ROOTDIR="${PWD}"
TOOLCHAIN_INSTALL_DIR="${TOOLCHAIN_ROOTDIR}/install"
TOOLCHAIN_SCRIPTS_DIR="${TOOLCHAIN_ROOTDIR}/scripts"
source "${TOOLCHAIN_SCRIPTS_DIR}/tool_kit.sh"

# ====================== Parameter parsing ======================
CLEAN_BUILD="__FALSE__"
BUILD_JOBS="$(get_nprocs)"
BUILD_SHARED_LIBS="ON"
DRY_RUN="__FALSE__"

show_help() {
  cat << EOF
Usage: $(basename "$0") [OPTIONS]

Generate CMake options from the toolchain configuration and build CP2K.

Options:
  -h, --help            Show this help message and exit
  -j N, -jN             Number of parallel build jobs
  --dry-run             Show generated CMake options only and then exit
  --clean               Remove the build directory before configuring, which means
                        rebuilding CP2K entirely
  --build-static        Set -DBUILD_SHARED_LIBS=OFF (default is ON)
EOF
}

while [ $# -ge 1 ]; do
  case "$1" in
    --dry-run)
      DRY_RUN="__TRUE__"
      ;;
    --clean)
      CLEAN_BUILD="__TRUE__"
      ;;
    -j)
      BUILD_JOBS="$2"
      shift
      ;;
    -j[0-9]*)
      BUILD_JOBS="${1#-j}"
      ;;
    --build-static)
      BUILD_SHARED_LIBS="OFF"
      ;;
    -h | --help)
      show_help
      exit 0
      ;;
    *)
      echo "ERROR: Unknown option: $1"
      show_help
      exit 1
      ;;
  esac
  shift
done

# ====================== Pre-checks ======================
if [ ! -f "${TOOLCHAIN_INSTALL_DIR}/setup" ]; then
  echo "Error: Toolchain is not installed. Please run ./install_cp2k_toolchain.sh first."
  exit 1
fi

# Load toolchain environment (required for with_xxx variables)
source "${TOOLCHAIN_INSTALL_DIR}/setup"
source "${TOOLCHAIN_INSTALL_DIR}/toolchain.conf"

printf "========================== %s =========================\n" \
  "Generating CMake options from toolchain"

# Require complete source tree
CP2K_ROOT=$(cd "${TOOLCHAIN_ROOTDIR}/../.." && pwd)
if [ -d "${CP2K_ROOT}/src" ]; then
  echo "Root directory of CP2K with source code is found as ${CP2K_ROOT}"
  echo "(path is exported to variable \${CP2K_ROOT})."
else
  report_error ${LINENO} "\${CP2K_ROOT} does not have subdirectory src."
  exit 1
fi
if [ -f "${CP2K_ROOT}/CMakeLists.txt" ] && [ -r "${CP2K_ROOT}/CMakeLists.txt" ]; then
  echo "\${CP2K_ROOT}/CMakeLists.txt exists; will be parsed for CMake options."
else
  report_error ${LINENO} "\${CP2K_ROOT}/CMakeLists.txt cannot be found or read."
  exit 1
fi
if [ -d "${CP2K_ROOT}/data" ]; then
  echo "Data directory ${CP2K_ROOT}/data is found and set as CP2K_DATA_DIR."
else
  report_error ${LINENO} "Data directory \${CP2K_ROOT}/data cannot be found."
  exit 1
fi

# Generate cmake options for compiling cp2k
CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=${CP2K_ROOT}/install -DCP2K_DATA_DIR=${CP2K_ROOT}/data -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
if [ -n "$(grep -- "--install-all" "${TOOLCHAIN_INSTALL_DIR}/setup")" ]; then
  CMAKE_OPTIONS+=" -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF"
  for toolchain_option in $(grep -i "dontuse" "${TOOLCHAIN_INSTALL_DIR}/toolchain.conf" | grep -vi "gcc" | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ "$(eval echo "\$with_${toolchain_option}")" != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' "${CP2K_ROOT}/CMakeLists.txt" | grep -i "${toolchain_option}" | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=OFF" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS+=" -D${ADDED_CMAKE_OPTION}=OFF"
      fi
    fi
  done
else
  # If MPI is used, set "CP2K_USE_MPI" to "ON"
  if [ "${MPI_MODE}" != "no" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_MPI=ON"
    if [ -n "$(grep "MPI_F08" "${TOOLCHAIN_INSTALL_DIR}/toolchain.env")" ]; then
      CMAKE_OPTIONS+=" -DCP2K_USE_MPI_F08=ON"
    fi
  fi
  # If GPU acceleration is used, add the option about GPU acceleration
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=${GPUVER}"
  elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=${GPUVER}"
  elif [ "${ENABLE_OPENCL}" = "__TRUE__" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_ACCEL=OPENCL"
  fi
  # Some options that should be specially considered:
  # Intel MKL includes FFTW
  if [ "${with_fftw}" != "__DONTUSE__" ] || [ "${MATH_MODE}" = "mkl" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_FFTW3=ON"
  fi
  # Mimic-MCL (MiMiC Communication Library)
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_MIMIC=ON"
  fi
  # Detect if any other dependencies is used and add the proper cmake option
  # Since "pugixml" and "gsl" are not mentioned in CMakeLists.txt, they will not be considered.
  for toolchain_option in $(grep -vi "dry\|list\|dontuse\|gcc\|cmake\|fftw\|mkl\|dbcsr" "${TOOLCHAIN_INSTALL_DIR}/toolchain.conf" | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ "$(eval echo "\$with_${toolchain_option}")" != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' "${CP2K_ROOT}/CMakeLists.txt" | grep -i "${toolchain_option}" | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=ON" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS+=" -D${ADDED_CMAKE_OPTION}=ON"
      fi
    fi
  done
fi

# Set build directory
BUILD_DIR="${CP2K_ROOT}/build"

# Show CMake options
echo "Generated CMake flags:"
for flag in ${CMAKE_OPTIONS}; do
  echo "   ${flag}"
done

if [ "${DRY_RUN}" != "__TRUE__" ]; then
  # ====================== Optional clean ======================
  if [ "${CLEAN_BUILD}" = "__TRUE__" ] && [ -d "${BUILD_DIR}" ]; then
    echo "Removing existing build directory: ${BUILD_DIR}"
    rm -rf "${BUILD_DIR}"
  fi

  mkdir -p "${BUILD_DIR}"

  # ====================== Configure ======================
  echo "================== CMake configuration ==================="
  echo "Source dir : ${CP2K_ROOT}"
  echo "Build  dir : ${BUILD_DIR}"
  echo "Shared libs: ${BUILD_SHARED_LIBS}"

  echo -n "Configuring ... "
  cmake -S "${CP2K_ROOT}" -B "${BUILD_DIR}" ${CMAKE_OPTIONS} > cmake.log 2>&1 || tail_excerpt cmake.log
  echo "done."
  echo "=========================================================="

  # ====================== Build ======================
  echo "==================== Building CP2K ======================="
  echo "Parallel jobs: ${BUILD_JOBS}"

  echo -n "Building CP2K ... "
  cmake --build "${BUILD_DIR}" --target install -j "${BUILD_JOBS}" > build.log 2>&1 || tail_excerpt build.log
  echo "done."
  echo "=========================================================="

  # Export variable for CMake options to cp2k_env file
  cat << EOF > "${CP2K_ROOT}/install/cp2k_env"
#!/bin/bash
export CP2K_ROOT="${CP2K_ROOT}"
source ${TOOLCHAIN_INSTALL_DIR}/setup
prepend_path PATH "${CP2K_ROOT}/install/bin"
prepend_path LD_LIBRARY_PATH "${CP2K_ROOT}/install/lib"
prepend_path PKG_CONFIG_PATH "${CP2K_ROOT}/install/lib/pkgconfig"
EOF

  echo
  echo "Done! Installed binaries are now available in: ${CP2K_ROOT}/install/bin"
  echo
  echo "Please always source the below script to load CP2K environment before running CP2K:"
  echo "  source ${CP2K_ROOT}/install/cp2k_env"
  echo
  echo "To run regtests:"
  echo "  ${CP2K_ROOT}/tests/do_regtest.py ${CP2K_ROOT}/install/bin psmp"
  echo "Run \`${CP2K_ROOT}/tests/do_regtest.py --help\` for available options of regtesting."
else
  echo "Since you run this script with \"--dry-run\", it now exits."
  echo "To build CP2K, drop off this flag."
  exit 0
fi

#EOF
