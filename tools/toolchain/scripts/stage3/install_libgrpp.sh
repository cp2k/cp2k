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

[ -f "${BUILDDIR}/setup_libgrpp" ] && rm "${BUILDDIR}/setup_libgrpp"

if [ "$with_libgrpp" != "__DONTUSE__" ]; then
  echo "==================== Using libgrpp ===================="

  cat << EOF >> "${BUILDDIR}/setup_libgrpp"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBGRPP"
EOF
fi

load "${BUILDDIR}/setup_libgrpp"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libgrpp"
