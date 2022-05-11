#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

export MATH_CFLAGS=''
export MATH_LDFLAGS=''
export MATH_LIBS=''

write_toolchain_env "${INSTALLDIR}"

case "$MATH_MODE" in
  mkl)
    "${SCRIPTDIR}"/stage2/install_mkl.sh "${with_mkl}"
    load "${BUILDDIR}/setup_mkl"
    ;;
  acml)
    "${SCRIPTDIR}"/stage2/install_acml.sh "${with_acml}"
    load "${BUILDDIR}/setup_acml"
    ;;
  openblas)
    "${SCRIPTDIR}"/stage2/install_openblas.sh "${with_openblas}"
    load "${BUILDDIR}/setup_openblas"
    ;;
  cray)
    # note the space is intentional so that the variable is
    # non-empty and can pass require_env checks
    export MATH_LDFLAGS="${MATH_LDFLAGS} "
    export MATH_LIBS="${MATH_LIBS} ${CRAY_EXTRA_LIBS}"
    ;;
esac

export CP_CFLAGS="${CP_CFLAGS} ${MATH_CFLAGS}"
export CP_LDFLAGS="${CP_LDFLAGS} ${MATH_LDFLAGS}"
export CP_LIBS="${CP_LIBS} ${MATH_LIBS}"

write_toolchain_env "${INSTALLDIR}"

#EOF
