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

export CP2K_ROOT=$(cd ${ROOTDIR}/../.. && pwd)

# ------------------------------------------------------------------------
# generate cmake options for compiling cp2k
# ------------------------------------------------------------------------
# Build the program in source tree for convenience
CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=../install -DCP2K_DATA_DIR=${CP2K_ROOT}/data"
if [ -n "$(grep -- "--install-all" ${INSTALLDIR}/setup)" ]; then
  CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF"
  # Since "--install-all" can be used together with "--with-PKG=no", an extra safeguard is added here
  for toolchain_option in $(grep -i "dontuse" ${INSTALLDIR}/toolchain.conf | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/target_compile_definitions(/,/)/p' ${CP2K_ROOT}/src/CMakeLists.txt | grep -i $toolchain_option | cut -d'{' -f2 | cut -d'}' -f1 | head -n 1)
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
    if [ -n $(grep "IF_MPI(-D__MPI_F08|)" "${INSTALLDIR}"/toolchain.conf) ]; then
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
  # When user chooses dftd4 rather than tblite in toolchain, the ADDED_CMAKE_OPTION below is set to CP2K_USE_TBLITE by mistake.
  if [ "${with_dftd4}" != "__DONTUSE__" ] && [ "${with_tblite}" = "__DONTUSE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_DFTD4=ON"
  fi
  # There is only MCL (MiMiC Communication Library) and no MiMiC that can be installed via toolchain, but if one chooses to install MCL then MiMiC should have been installed?
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MIMIC=ON"
  fi
  # Detect if any other dependencies is used and add the proper cmake option. Since "pugixml" and "gsl" are not mentioned in CMakeLists.txt, they will not be considered.
  for toolchain_option in $(grep -vi "dontuse\|list\|cmake\|fftw\|dbcsr\|dftd4" ${INSTALLDIR}/toolchain.conf | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/target_compile_definitions(/,/)/p' ${CP2K_ROOT}/src/CMakeLists.txt | grep -i $toolchain_option | cut -d'{' -f2 | cut -d'}' -f1 | head -n 1)
      # Use "if-then" below can avoid generating empty "-D=ON" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -D${ADDED_CMAKE_OPTION}=ON"
      fi
    fi
  done
fi

# -------------------------
# print out user instructions
# -------------------------

cat << EOF

========================== usage =========================
Done! The "build" directory can now be removed.

To use the installed tools and libraries and cp2k version compiled with it you will first need to execute at the prompt:
  source ${SETUPFILE}

It's recommended for you to build CP2K like this after executing above command:
  cd ${CP2K_ROOT}
  mkdir build && cd build
  cmake .. ${CMAKE_OPTIONS}
  make install -j $(get_nprocs)

When completed, you can run "make clean" or delete this build directory to free up some space.
EOF

#EOF
