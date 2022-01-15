#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

for ii in $tool_list; do
  load "${BUILDDIR}/setup_${ii}"
done

# ------------------------------------------------------------------------
# Install or compile packages using newly installed tools
# ------------------------------------------------------------------------

# setup compiler flags, leading to nice stack traces on crashes but still optimised
if [ "${generic}" = "__TRUE__" ]; then
  CFLAGS="-O2 -fPIC -fno-omit-frame-pointer -g -mtune=generic ${TSANFLAGS}"
else
  CFLAGS="-O2 -fPIC -fno-omit-frame-pointer -g -march=native -mtune=native ${TSANFLAGS}"
fi
FFLAGS="${CFLAGS} -fbacktrace"
CXXFLAGS="${CFLAGS}"
F77FLAGS="${FFLAGS}"
F90FLAGS="${FFLAGS}"
FCFLAGS="${FFLAGS}"

export CFLAGS=$(allowed_gcc_flags ${CFLAGS})
export FFLAGS=$(allowed_gfortran_flags ${FFLAGS})
export F77FLAGS=$(allowed_gfortran_flags ${F77FLAGS})
export F90FLAGS=$(allowed_gfortran_flags ${F90FLAGS})
export FCFLAGS=$(allowed_gfortran_flags ${FCFLAGS})
export CXXFLAGS=$(allowed_gxx_flags ${CXXFLAGS})
export LDFLAGS="${TSANFLAGS}"

# get system arch information using OpenBLAS prebuild
"${SCRIPTDIR}"/get_openblas_arch.sh
load "${BUILDDIR}/openblas_arch"

write_toolchain_env "${INSTALLDIR}"

#EOF
