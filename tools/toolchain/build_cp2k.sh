#!/bin/bash -e

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

TOOLCHAIN_ROOTDIR="${PWD}"
# Exit this script if it is not called from the ./tools/toolchain directory
if [ "${TOOLCHAIN_ROOTDIR}" != "${SCRIPT_DIR}" ]; then
  cat << EOF
ERROR: Incorrect execution location.
The absolute path of the build_cp2k.sh script is at:
  ${SCRIPT_DIR}
Actual working directory where it is currently called:
  ${TOOLCHAIN_ROOTDIR}
Please enter the absolute path above before executing the build_cp2k.sh script
so that subsequent scripts can be found and files can be placed correctly.
EOF
  exit 1
fi
TOOLCHAIN_SCRIPTS_DIR="${TOOLCHAIN_ROOTDIR}/scripts"
source "${TOOLCHAIN_ROOTDIR}/toolchain_settings"
source "${TOOLCHAIN_SCRIPTS_DIR}/tool_kit.sh"

# ====================== Parameter parsing ======================
CP2K_ROOT=$(cd "${TOOLCHAIN_ROOTDIR}/../.." && pwd)
CMAKE_INSTALL_PREFIX=${CP2K_ROOT}/install
CLEAN_BUILD="__FALSE__"
DEBUG_BUILD="__FALSE__"
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
  --prefix              Set CMAKE_INSTALL_PREFIX (default is ${CP2K_ROOT}/install)
  --dry-run             Show generated CMake options only and then exit
  --clean               Remove the build directory before configuring, which means
                        rebuilding CP2K entirely
  --debug               Build debug version of CP2K (-DCMAKE_BUILD_TYPE=Debug)
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
    --prefix)
      if [[ "${2}" != /* ]]; then
        report_error "The path for --prefix must be an absolute path."
        exit 1
      fi
      CMAKE_INSTALL_PREFIX="${2}"
      shift
      ;;
    --debug)
      DEBUG_BUILD="__TRUE__"
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
# Require complete source tree
# A minimum working environment for this script should be as follows (assuming out-of-tree build not required):
# cp2k                                  <- variable ${CP2K_ROOT}; CMake option -S
# ├── CMakeLists.txt                    <- file to be parsed for generating cmake options
# ├── cmake                             <- directory containing CMake files
# ├── data                              <- CP2K data directory; CMake option -DCP2K_DATA_DIR\
# ├── build                             <- to-be-created; CMake option -B
# ├── install                           <- to-be-created; CMake option -DCMAKE_INSTALL_PREFIX
# ├── src                               <- CP2K source code directory
# └── tools
#     └── toolchain                     <- working directory; variable ${TOOLCHAIN_ROOTDIR}
#         ├── build_cp2k.sh             <- this script
#         ├── install_cp2k_toolchain.sh <- script being executed before calling this script
#         ├── scripts                   <- directory with toolchain scripts; ${TOOLCHAIN_SCRIPTS_DIR}
#         │   └── tool_kit.sh
#         └── install                   <- directory with installed dependencies; ${TOOLCHAIN_INSTALL_DIR}
#             ├── setup                 <- * file to be used for building CP2K
#             ├── toolchain.conf        <- * file to be parsed for generating cmake options
#             └── toolchain.env         <- * file to be parsed for generating cmake options (MPI_F08)
#
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
  echo "Data directory ${CP2K_ROOT}/data is found."
else
  report_error ${LINENO} "Data directory \${CP2K_ROOT}/data cannot be found."
  exit 1
fi

# Require finished toolchain
if [ ! -f "${TOOLCHAIN_INSTALL_DIR}/setup" ]; then
  echo "Error: Toolchain is not installed. Please run ./install_cp2k_toolchain.sh first."
  exit 1
fi

# Disallow combination of installing toolchain outside the source tree and CP2K under the source tree
if [ "${TOOLCHAIN_INSTALL_DIR}" != "${TOOLCHAIN_ROOTDIR}"/install ] &&
  [[ ${CMAKE_INSTALL_PREFIX} == ${CP2K_ROOT}/* ]]; then
  echo
  cat << EOF
ERROR: You toolchain installation is outside the source tree but the install
prefix of CP2K is under the source tree, which is disallowed by this script.
The script will now abort; please manually set "--prefix" to a proper path.
EOF
  exit 1
fi

# Load toolchain environment (required for with_xxx variables)
source "${TOOLCHAIN_INSTALL_DIR}/setup"
source "${TOOLCHAIN_INSTALL_DIR}/toolchain.conf"

# Generate cmake options for compiling cp2k
CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}"
if [[ ${CMAKE_INSTALL_PREFIX} == ${CP2K_ROOT}/* ]]; then
  CMAKE_OPTIONS+=" -DCP2K_DATA_DIR=${CP2K_ROOT}/data"
fi
if [ ${DEBUG_BUILD} == "__TRUE__" ]; then
  CMAKE_OPTIONS+=" -DCMAKE_BUILD_TYPE=Debug"
fi
if [ -n "$(grep -- "--install-all" "${TOOLCHAIN_INSTALL_DIR}/setup")" ]; then
  CMAKE_OPTIONS+=" -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF"
  for toolchain_option in $(grep -i "dontuse" "${TOOLCHAIN_INSTALL_DIR}/toolchain.conf" |
    grep -Evi "gcc|amd|intel" | cut -d'_' -f2 | cut -d'=' -f1); do
    var_name="with_${toolchain_option}"
    if [ "${!var_name}" != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' "${CP2K_ROOT}/CMakeLists.txt" |
        grep -i "${toolchain_option}" | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=OFF" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS+=" -D${ADDED_CMAKE_OPTION}=OFF"
      fi
    fi
  done
else
  # If MPI is used, set "CP2K_USE_MPI" to "ON"
  if [ "${mpi_mode}" != "no" ]; then
    CMAKE_OPTIONS+=" -DCP2K_USE_MPI=ON -DCP2K_USE_MPI_F08=ON"
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
  for toolchain_option in $(grep "with" "${TOOLCHAIN_INSTALL_DIR}"/toolchain.conf |
    grep -Evi "dontuse|gcc|amd|intel|cmake|fftw|mkl|dbcsr" | cut -d'_' -f2 | cut -d'=' -f1); do
    var_name="with_${toolchain_option}"
    if [ "${!var_name}" != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' "${CP2K_ROOT}/CMakeLists.txt" |
        grep -i "${toolchain_option}" | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=ON" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS+=" -D${ADDED_CMAKE_OPTION}=ON"
      fi
    fi
  done
fi
# If GPU acceleration is used, add the option about GPU acceleration
if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
  CMAKE_OPTIONS+=" -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=${GPU_VER}"
elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
  CMAKE_OPTIONS+=" -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=${GPU_VER}"
elif [ "${ENABLE_OPENCL}" = "__TRUE__" ]; then
  CMAKE_OPTIONS+=" -DCP2K_USE_ACCEL=OPENCL"
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
  echo "Install dir: ${CMAKE_INSTALL_PREFIX}"
  echo "Shared libs: ${BUILD_SHARED_LIBS}"

  echo -n "Configuring ... "
  cmake -S "${CP2K_ROOT}" -B "${BUILD_DIR}" ${CMAKE_OPTIONS} > cmake.log 2>&1 || tail_excerpt cmake.log
  echo "done."

  # ====================== Build ======================
  echo "==================== Building CP2K ======================="
  echo "Parallel jobs: ${BUILD_JOBS}"

  echo -n "Building CP2K ... "
  if ! cmake --build "${BUILD_DIR}" --target install -j "${BUILD_JOBS}" > build.log 2>&1; then
    echo "failed."
    tail_excerpt build.log
    exit 1
  else
    echo "done."
  fi
  echo "=========================================================="

  # Export variable for CMake options to cp2k_env file
  cat << EOF > "${CMAKE_INSTALL_PREFIX}/cp2k_env"
#!/bin/bash
source ${TOOLCHAIN_INSTALL_DIR}/setup
prepend_path PATH "${CMAKE_INSTALL_PREFIX}/bin"
prepend_path LD_LIBRARY_PATH "${CMAKE_INSTALL_PREFIX}/lib"
prepend_path PKG_CONFIG_PATH "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}"
EOF
  cat << EOF
Done! Installed binaries are now available in: ${CMAKE_INSTALL_PREFIX}/bin

It's suggested to run regtests after installation:
  ${CP2K_ROOT}/tests/do_regtest.py ${CMAKE_INSTALL_PREFIX}/bin psmp
Run \`${CP2K_ROOT}/tests/do_regtest.py --help\` for help message.

Please always source this script to load CP2K environment before running CP2K:
  source ${CMAKE_INSTALL_PREFIX}/cp2k_env

If you want to clean the build cache (except cached CMake files) after
installation, run:
  cmake --build "${BUILD_DIR}" --target clean
EOF
else
  cat << EOF
Since you run this script with \"--dry-run\", it now exits.
To build CP2K, drop off this flag and re-run this script.
EOF
fi

#EOF
