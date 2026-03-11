#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

# This script assumes a working environment as follows:
# cp2k                                  <- variable ${CP2K_ROOT}
# ├── CMakeLists.txt                    <- * file to be parsed in this script
# ├── data                              <- CMake option -DCP2K_DATA_DIR
# ├── build                             <- to-be-created CP2K CMake build directory
# ├── install                           <- to-be-created CP2K install directory
# ├── src                               <- CP2K source code directory
# └── tools
#     └── toolchain                     <- working directory; variable ${ROOTDIR}
#         ├── install_cp2k_toolchain.sh <- script being executed calling this script
#         ├── scripts                   <- variable ${SCRIPT_DIR}
#         │   ├── common_vars.sh
#         │   ├── tool_kit.sh
#         │   └── generate_cmake_options.sh <- this script
#         └── install                   <- variable ${INSTALLDIR}
#             ├── setup                 <- * file to be parsed in this script
#             ├── toolchain.conf        <- * file to be parsed in this script
#             └── toolchain.env         <- * file to be parsed in this script
#
# First, validate existence of relevant upper-level directory and file
printf "\n========================== %s =========================\n" \
  "Generating CMake options for building CP2K"
export CP2K_ROOT=$(cd "${ROOTDIR}/../.." && pwd)
cat << EOF
Root directory of CP2K is assumed to be ${CP2K_ROOT} (variable \${CP2K_ROOT}).
Build directory will be ${CP2K_ROOT}/build for executing CMake.
Install directory will be ${CP2K_ROOT}/install for CP2K binaries, libraries etc.
EOF
CMAKE_OPTIONS="-S ${CP2K_ROOT} -B ${CP2K_ROOT}/build"
CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCMAKE_INSTALL_PREFIX=${CP2K_ROOT}/install"
CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCMAKE_INSTALL_LIBDIR=lib"
CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCMAKE_SKIP_RPATH=ON"
if [ -d "${CP2K_ROOT}/data" ]; then
  echo "Data directory ${CP2K_ROOT}/data is found and set as -DCP2K_DATA_DIR."
  CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_DATA_DIR=${CP2K_ROOT}/data"
else
  report_warning ${LINENO} "Data directory ${CP2K_ROOT}/data cannot be found;
please set -DCP2K_DATA_DIR to actual path manually."
fi
if [ -f "${CP2K_ROOT}/CMakeLists.txt" ] && [ -r "${CP2K_ROOT}/CMakeLists.txt" ]; then
  echo "${CP2K_ROOT}/CMakeLists.txt exists; will be parsed for CMake options."
else
  report_error ${LINENO} "${CP2K_ROOT}/CMakeLists.txt cannot be found or read;
suggested CMake option will be incomplete and/or incorrect."
  return 1
fi

# ------------------------------------------------------------------------
# generate cmake options for compiling cp2k
# ------------------------------------------------------------------------
# Build the program in source tree for convenience
if [ -n "$(grep -- "--install-all" ${INSTALLDIR}/setup)" ]; then
  CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF"
  # Since "--install-all" can be used together with "--with-PKG=no", an extra safeguard is added here
  for toolchain_option in $(grep -i "dontuse" ${INSTALLDIR}/toolchain.conf | grep -vi "gcc" | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' ${CP2K_ROOT}/CMakeLists.txt | grep -i $toolchain_option | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=OFF" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -D${ADDED_CMAKE_OPTION}=OFF"
      fi
    fi
  done
else
  # If MPI is used, set "CP2K_USE_MPI" to "ON"
  if [ "${MPI_MODE}" != "no" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MPI=ON"
    if [ -n "$(grep "MPI_F08" "${INSTALLDIR}"/toolchain.env)" ]; then
      CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MPI_F08=ON"
    fi
  fi
  # If GPU acceleration is used, add the option about GPU acceleration
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=${GPUVER}"
  elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=${GPUVER}"
  elif [ "${ENABLE_OPENCL}" = "__TRUE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=OPENCL"
  fi
  # Some options that should be specially considered:
  # Intel MKL includes FFTW
  if [ "${with_fftw}" != "__DONTUSE__" ] || [ "${MATH_MODE}" = "mkl" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_FFTW3=ON"
  fi
  # There is only MCL (MiMiC Communication Library) and no MiMiC that can be installed via toolchain, but if one chooses to install MCL then MiMiC should have been installed?
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MIMIC=ON"
  fi
  # Detect if any other dependencies is used and add the proper cmake option. Since "pugixml" and "gsl" are not mentioned in CMakeLists.txt, they will not be considered.
  for toolchain_option in $(grep -vi "dry\|list\|dontuse\|gcc\|cmake\|fftw\|mkl\|dbcsr" ${INSTALLDIR}/toolchain.conf | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/option(/,/)/p' ${CP2K_ROOT}/CMakeLists.txt | grep -i $toolchain_option | awk '{print $1}' | cut -d'(' -f2 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=ON" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -D${ADDED_CMAKE_OPTION}=ON"
      fi
    fi
  done
fi

# Export variable for CMake options to setup file
printf 'export CP2K_CMAKE_OPTIONS="%s"\n' "${CMAKE_OPTIONS}" >> "${SETUPFILE}"
cat << EOF

CMake options from parsing files are now collected in the variable
\${CP2K_CMAKE_OPTIONS}, which is exported at the end of setup file:
${SETUPFILE}
EOF

# -------------------------
# print out user instructions
# -------------------------
if [ "${dry_run}" = "__TRUE__" ]; then
  cat << EOF
Suggested cmake command if toolchain is built with your options:
  cmake ${CMAKE_OPTIONS}
EOF
else
  echo
  cat << EOF | tee ${INSTALLDIR}/cp2k_installation_guide.txt
========================== Epilogue =========================
Toolchain is now ready for building CP2K! Instructions for next steps:

(1) Optional - remove packages in ./build directory to free up disk space:
      rm -rf ${BUILDDIR}
    However, do NOT delete or move the ./install directory from now on.

(2) Required - source setup file to activate toolchain-configured dependencies:
      source ${SETUPFILE}
    This setup file MUST also be sourced whenever CP2K built with this toolchain
    is executed. If modules have been used to estabilish environment variables
    and paths, remember to load these modules prior to sourcing setup file.

(3) Recommended - go to root directory of CP2K and build executable:
      cd ${CP2K_ROOT}
      mkdir build && cd build
      cmake ${CMAKE_OPTIONS}
      make install -j $(get_nprocs)
    The CMake command shows the options from generate_cmake_options.sh. For more
    CMake options, see ${CP2K_ROOT}/CMakeLists.txt.
    Alternatively the CMake command may take the same options from a variable as
    exported in the setup file, no need to copy-paste long lines in terminal:
      cmake "\${CP2K_CMAKE_OPTIONS}"
    It may be helpful to also save a copy of command line messages to log files:
      cmake "\${CP2K_CMAKE_OPTIONS}" 2>&1 | tee make.log
      make install -j $(get_nprocs) 2>&1 | tee install.log

(4) Optional - once build is completed, remove ./build directory like in (1):
      make clean
    Again, do not move the ./install directory from now on.

For more information about available build options, see:
https://manual.cp2k.org/trunk/getting-started/build-from-source.html.
This message is saved to "${INSTALLDIR}/cp2k_installation_guide.txt".
EOF
fi

#EOF
