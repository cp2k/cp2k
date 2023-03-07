#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

source ${SCRIPT_DIR}/common_vars.sh
source ${SCRIPT_DIR}/tool_kit.sh
source ${SCRIPT_DIR}/signal_trap.sh
source ${INSTALLDIR}/toolchain.conf
source ${INSTALLDIR}/toolchain.env

for ii in $tool_list; do
  load "${BUILDDIR}/setup_${ii}"
done

# ------------------------------------------------------------------------
# Install or compile packages using newly installed tools
# ------------------------------------------------------------------------

# Setup compiler flags, leading to nice stack traces on crashes but still optimised
if [ "${with_intel}" != "__DONTUSE__" ]; then
  CFLAGS="-O2 -fPIC -fp-model=precise -funroll-loops -g -qopenmp -qopenmp-simd -traceback"
  if [ "${TARGET_CPU}" = "native" ]; then
    CFLAGS="${CFLAGS} -xHost"
  elif [ "${TARGET_CPU}" = "generic" ]; then
    CFLAGS="${CFLAGS} -mtune=${TARGET_CPU}"
  else
    CFLAGS="${CFLAGS} -march=${TARGET_CPU} -mtune=${TARGET_CPU}"
  fi
  FFLAGS="${CFLAGS}"
else
  CFLAGS="-O2 -fPIC -fno-omit-frame-pointer -fopenmp -g"
  if [ "${TARGET_CPU}" = "generic" ]; then
    CFLAGS="${CFLAGS} -mtune=generic ${TSANFLAGS}"
  else
    CFLAGS="${CFLAGS} -march=${TARGET_CPU} -mtune=${TARGET_CPU} ${TSANFLAGS}"
  fi
  FFLAGS="${CFLAGS} -fbacktrace"
fi
CXXFLAGS="${CFLAGS}"
F77FLAGS="${FFLAGS}"
F90FLAGS="${FFLAGS}"
FCFLAGS="${FFLAGS}"

if [ "${with_intel}" == "__DONTUSE__" ]; then
  export CFLAGS="$(allowed_gcc_flags ${CFLAGS})"
  export FFLAGS="$(allowed_gfortran_flags ${FFLAGS})"
  export F77FLAGS="$(allowed_gfortran_flags ${F77FLAGS})"
  export F90FLAGS="$(allowed_gfortran_flags ${F90FLAGS})"
  export FCFLAGS="$(allowed_gfortran_flags ${FCFLAGS})"
  export CXXFLAGS="$(allowed_gxx_flags ${CXXFLAGS})"
else
  # TODO Check functions for allowed Intel compiler flags
  export CFLAGS
  export FFLAGS
  export F77FLAGS
  export F90FLAGS
  export FCFLAGS
  export CXXFLAGS
fi
export LDFLAGS="${TSANFLAGS}"

# get system arch information using OpenBLAS prebuild
${SCRIPTDIR}/get_openblas_arch.sh
load "${BUILDDIR}/openblas_arch"

write_toolchain_env "${INSTALLDIR}"

#EOF
