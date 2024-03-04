#!/bin/bash -e

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libsmeagol_ver="1.2"
libsmeagol_sha256="aaf7f2238182413d55b6ad819e02b90c95db7b9f7f9c22433085ce73b5f835bd"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libsmeagol" ] && rm "${BUILDDIR}/setup_libsmeagol"

LIBSMEAGOL_CFLAGS=""
LIBSMEAGOL_LDFLAGS=""
LIBSMEAGOL_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libsmeagol" in
  __INSTALL__)
    echo "==================== Installing SMEAGOL ===================="
    echo "__INSTALL__ is not supported; please manually install libsmeagol"
    echo "and rerun the toolchain using --with-libsmeagol=<path> option"
    exit 1
    ;;

  __SYSTEM__)
    echo "==================== Finding SMEAGOL from system paths ===================="
    echo "Please rerun the toolchain using --with-libsmeagol=<path> option"
    exit 1
    ;;

  __DONTUSE__) ;;

  *)
    echo "==================== Linking SMEAGOL to user paths ===================="
    pkg_install_dir="$with_libsmeagol"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/obj"
    LIBSMEAGOL_CFLAGS="-I'${pkg_install_dir}/obj'"
    LIBSMEAGOL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libsmeagol" != "__DONTUSE__" ]; then
  LIBSMEAGOL_LIBS="-lsmeagol"
  if [ "$with_libsmeagol" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libsmeagol"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
    cat "${BUILDDIR}/setup_libsmeagol" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_libsmeagol"
export LIBSMEAGOL_CFLAGS="${LIBSMEAGOL_CFLAGS}"
export LIBSMEAGOL_LDFLAGS="${LIBSMEAGOL_LDFLAGS}"
export LIBSMEAGOL_LIBS="${LIBSMEAGOL_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__SMEAGOL"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBSMEAGOL_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBSMEAGOL_LDFLAGS}"
export CP_LIBS="${LIBSMEAGOL_LIBS} \${CP_LIBS}"
export LIBSMEAGOL_ROOT="$pkg_install_dir"
EOF
fi

load "${BUILDDIR}/setup_libsmeagol"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libsmeagol"
