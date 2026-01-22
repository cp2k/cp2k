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
  # If MPI used, set "CP2K_USE_MPI" to "ON"
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
  # Detect if any other dependencies is used simply via the variables in file "${INSTALLDIR}"/toolchain.env
  # Disadvantage: If a new dependencies was added to toolchain, the option will have to be manually added to this script.
  # Example:
  #  if [ "${with_abc}" != "__DONTUSE__" ]; then
  #    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ABC=ON"
  #  fi
  if [ "${with_fftw}" != "__DONTUSE__" ] || [ "MATH_MODE" = "mkl" ]; then
    # Intel MKL includes FFTW
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_FFTW3=ON"
  fi
  if [ "${with_libint}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_LIBINT2=ON"
  fi
  if [ "${with_libxc}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_LIBXC=ON"
  fi
  if [ "${with_libxsmm}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_LIBXSMM=ON"
  fi
  if [ "${with_sirius}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_SIRIUS=ON"
  fi
  if [ "${with_spglib}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_SPGLIB=ON"
  fi
  if [ "${with_cosma}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_COSMA=ON"
  fi
  if [ "${with_elpa}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ELPA=ON"
  fi
  if [ "${with_hdf5}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_HDF5=ON"
  fi
  if [ "${with_libvori}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_VORI=ON"
  fi
  if [ "${with_libvdwxc}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_LIBVDWXC=ON"
  fi
  if [ "${with_dftd4}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_DFTD4=ON"
  fi
  if [ "${with_tblite}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_TBLITE=ON"
    if [ -f ${INSTALLDIR}/tblite-${TBLITE_VER}/bin/dftd4 ] && [ -d ${INSTALLDIR}/tblite-${TBLITE_VER}/include/dftd4 ]; then
      export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_DFTD4=ON"
    fi
  fi
  if [ "${with_plumed}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_PLUMED=ON"
  fi
  if [ "${with_spla}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_SPLA=ON"
  fi
  if [ "${with_trexio}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_TREXIO=ON"
  fi
  if [ "${with_cusolvermp}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_CUSOLVER_MP=ON"
  fi
  if [ "${with_libtorch}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_LIBTORCH=ON"
  fi
  if [ "${with_deepmd}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_DEEPMD=ON"
  fi
  if [ "${with_libsmeagol}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_LIBSMEAGOL=ON"
  fi
  if [ "${with_greenx}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_GREENX=ON"
  fi
  if [ "${with_ace}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_ACE=ON"
  fi
  if [ "${with_mcl}" != "__DONTUSE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_MIMIC=ON"
  fi
  if [ "${enable_cray}" = "__TRUE__" ]; then
    export CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCP2K_USE_CRAY_PM_ENERGY=ON"
  fi
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
