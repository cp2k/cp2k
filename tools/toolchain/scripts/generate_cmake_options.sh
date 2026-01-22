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
export CMAKE_OPTIONS="-DCMAKE_INSTALL_PREFIX=../install"
export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_DATA_DIR=${CP2K_ROOT}/data"
if [ -n "$(grep -- "--install-all" ${INSTALLDIR}/setup)" ]; then
  export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_EVERYTHING=ON -DCP2K_USE_DLAF=OFF -DCP2K_USE_PEXSI=OFF"
else
  # If MPI is used, set "CP2K_USE_MPI" to "ON"
  if [ "${MPI_MODE}" != "no" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MPI=ON"
  fi
  # If GPU acceleration is used, add the option about GPU acceleration
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=CUDA -DCP2K_WITH_GPU=${GPUVER}"
  else
    if [ "${ENABLE_HIP}" = "__TRUE__" ]; then
      export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=HIP -DCP2K_WITH_GPU=${GPUVER}"
    else
      if [ "${ENABLE_OPENCL}" = "__TRUE__" ]; then
        export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACCEL=OPENCL"
      fi
    fi
  fi
  # Some options that should be specially considered:
  # Intel MKL includes FFTW
  if [ "${with_fftw}" != "__DONTUSE__" ] || [ "MATH_MODE" = "mkl" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_FFTW3=ON"
  fi
  # There is only MCL (MiMiC Communication Library) and no MiMiC that can be installed via toolchain, but if one chooses to install MCL then MiMiC should have been installed?
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MIMIC=ON"
  fi
  # Detect if any other dependencies is used and add the proper cmake option. Since "pugixml" and "gsl" are not mentioned in CMakeLists.txt, they will not be considered.
  for toolchain_option in $(grep -vi "dontuse\|list\|cmake\|fftw\|dbcsr" ${INSTALLDIR}/toolchain.conf | cut -d'_' -f2 | cut -d'=' -f1); do
    if [ $(eval echo "$with_"'$toolchain_option') != "__DONTUSE__" ]; then
      ADDED_CMAKE_OPTION=$(sed -n '/target_compile_definitions(/,/)/p' ${CP2K_ROOT}/src/CMakeLists.txt | grep -i $toolchain_option | cut -d'{' -f2 | cut -d'}' -f1 | grep -i $toolchain_option | head -n 1)
      # The reason that "grep -i $toolchain_option" appears twice is related to DFT-D4 which is a submodule of tblite: if user choose dftd4 rather than tblite in toolchain, the option will be added as CP2K_USE_TBLITE by mistake without the second "grep" command.
      # Use "if-then" below can avoid generating empty "-D=ON" options
      if [ -n "${ADDED_CMAKE_OPTION}" ]; then
        export CMAKE_OPTIONS="${CMAKE_OPTIONS} -D${ADDED_CMAKE_OPTION}=ON"
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
