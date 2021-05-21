#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

# math core libraries, need to use reflapack for valgrind builds, as
# many fast libraries are not necessarily valgrind clean
export REF_MATH_CFLAGS=''
export REF_MATH_LDFLAGS=''
export REF_MATH_LIBS=''
export FAST_MATH_CFLAGS=''
export FAST_MATH_LDFLAGS=''
export FAST_MATH_LIBS=''

write_toolchain_env "${INSTALLDIR}"

"${SCRIPTDIR}"/stage2/install_reflapack.sh "${with_reflapack}"
load "${BUILDDIR}/setup_reflapack"

case "$FAST_MATH_MODE" in
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
    export FAST_MATH_LDFLAGS="${FAST_MATH_LDFLAGS} "
    export FAST_MATH_LIBS="${FAST_MATH_LIBS} ${CRAY_EXTRA_LIBS}"
    ;;
esac

if [ "$with_reflapack" = "__DONTUSE__" ]; then
  # if we don't build the refereence blas/lapack implementations,
  # make sure we still link against a BLAS/LAPACK implementation in *dbg profiles
  export REF_MATH_CFLAGS="${FAST_MATH_CFLAGS}"
  export REF_MATH_LDFLAGS="${FAST_MATH_LDFLAGS}"
  export REF_MATH_LIBS="${FAST_MATH_LIBS}"
fi

export MATH_CFLAGS="${FAST_MATH_CFLAGS}"
export MATH_LDFLAGS="${FAST_MATH_LDFLAGS}"
export MATH_LIBS="${FAST_MATH_LIBS}"

export CP_CFLAGS="${CP_CFLAGS} IF_DEBUG(${REF_MATH_CFLAGS}|IF_VALGRIND(${REF_MATH_CFLAGS}|${FAST_MATH_CFLAGS}))"
export CP_LDFLAGS="${CP_LDFLAGS} IF_DEBUG(${REF_MATH_LDFLAGS}|IF_VALGRIND(${REF_MATH_LDFLAGS}|${FAST_MATH_LDFLAGS}))"
export CP_LIBS="${CP_LIBS} IF_DEBUG(${REF_MATH_LIBS}|IF_VALGRIND(${REF_MATH_LIBS}|${FAST_MATH_LIBS}))"

write_toolchain_env "${INSTALLDIR}"

#EOF
