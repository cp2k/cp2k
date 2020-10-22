#!/bin/bash -e

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

for ii in $tool_list ; do
    load "${BUILDDIR}/setup_${ii}"
done

# ------------------------------------------------------------------------
# Install or compile packages using newly installed tools
# ------------------------------------------------------------------------

# setup compiler flags, leading to nice stack traces on crashes but
# still optimised
CFLAGS="-O2 -fno-omit-frame-pointer -g -march=native -mtune=native ${TSANFLAGS}"
CXXFLAGS="${CFLAGS}"
FFLAGS="-O2 -fbacktrace -fno-omit-frame-pointer -g -march=native -mtune=native ${TSANFLAGS}"
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
"${SCRIPTDIR}"/get_openblas_arch.sh; load "${BUILDDIR}/openblas_arch"

write_toolchain_env "${INSTALLDIR}"

#EOF
