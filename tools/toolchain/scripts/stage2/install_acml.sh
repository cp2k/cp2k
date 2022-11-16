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

[ -f "${BUILDDIR}/setup_acml" ] && rm "${BUILDDIR}/setup_acml"

ACML_CFLAGS=''
ACML_LDFLAGS=''
ACML_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_acml" in
  __INSTALL__)
    echo "==================== Installing ACML ===================="
    report_error $LINENO "To install ACML you should either contact your system administrator or go to https://developer.amd.com/tools-and-sdks/archive/amd-core-math-library-acml/acml-downloads-resources/ and download the correct version for your system."
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding ACML from system paths ===================="
    check_lib -lacml "ACML"
    add_include_from_paths ACML_CFLAGS "acml.h" $INCLUDE_PATHS
    add_lib_from_paths ACML_LDFLAGS "libacml.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking ACML to user paths ===================="
    pkg_install_dir="$with_acml"
    check_dir "${pkg_install_dir}/lib"
    ACML_CFLAGS="-I'${pkg_install_dir}/include'"
    ACML_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_acml" != "__DONTUSE__" ]; then
  ACML_LIBS="-lacml"
  if [ "$with_acml" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_acml"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_acml" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_acml"
export ACML_CFLAGS="${ACML_CFLAGS}"
export ACML_LDFLAGS="${ACML_LDFLAGS}"
export ACML_LIBS="${ACML_LIBS}"
export MATH_CFLAGS="\${MATH_CFLAGS} ${ACML_CFLAGS}"
export MATH_LDFLAGS="\${MATH_LDFLAGS} ${ACML_LDFLAGS}"
export MATH_LIBS="\${MATH_LIBS} ${ACML_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_acml"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "acml"
