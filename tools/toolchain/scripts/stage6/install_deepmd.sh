#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

deepmd_ver="2.2.7"
deepmd_pkg="libdeepmd_c-${deepmd_ver}.tar.gz"
deepmd_sha256="605bbb0c3bcc847ecbfe7326bac9ff4cac5690accc8083122e3735290bb923ae"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_deepmd" ] && rm "${BUILDDIR}/setup_deepmd"

DEEPMD_LDFLAGS=''
DEEPMD_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_deepmd" in
  __INSTALL__)
    echo "==================== Installing DeePMD ===================="
    pkg_install_dir="${INSTALLDIR}/libdeepmd_c-${deepmd_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    deepmd_root="${pkg_install_dir}"
    if verify_checksums "${install_lock_file}"; then
      echo "libdeepmd_c-${deepmd_ver} is already installed, skipping it."
    else
      if [ -f ${deepmd_pkg} ]; then
        echo "${deepmd_pkg} is found"
      else
        download_pkg_from_cp2k_org "${deepmd_sha256}" "${deepmd_pkg}"
      fi
      [ -d libdeepmd_c ] && rm -rf libdeepmd_c
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf ${deepmd_pkg}
      mv libdeepmd_c ${pkg_install_dir}
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    DEEPMD_DFLAGS="-D__DEEPMD"
    DEEPMD_CFLAGS="-I'${pkg_install_dir}/include'"
    DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_c -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding DeePMD from system paths ===================="
    check_lib -ldeepmd "DEEPMD"
    add_lib_from_paths DEEPMD_LDFLAGS "libdeepmd*" $LIB_PATHS
    add_include_from_paths DEEPMD_CFLAGS "deepmd" $INCLUDE_PATHS
    DEEPMD_DFLAGS="-D__DEEPMD"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking DEEPMD to user paths ===================="
    pkg_install_dir="$with_deepmd"
    check_dir "${pkg_install_dir}/include/deepmd"
    check_dir "${pkg_install_dir}/lib"
    DEEPMD_DFLAGS="-D__DEEPMD"
    DEEPMD_CFLAGS="-I'${pkg_install_dir}/include'"
    DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_c -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_deepmd" != "__DONTUSE__" ]; then
  DEEPMD_LIBS='-ldeepmd_c -lstdc++'
  if [ "$with_deepmd" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_deepmd"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
EOF
    cat "${BUILDDIR}/setup_deepmd" >> $SETUPFILE
  fi

  cat << EOF >> "${BUILDDIR}/setup_deepmd"
export DEEPMD_DFLAGS="${DEEPMD_DFLAGS}"
export DEEPMD_CFLAGS="${DEEPMD_CFLAGS}"
export DEEPMD_LDFLAGS="${DEEPMD_LDFLAGS}"
export DEEPMD_LIBS="${DEEPMD_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} ${DEEPMD_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${DEEPMD_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${DEEPMD_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${DEEPMD_LIBS}"
EOF
  cat << EOF >> "${INSTALLDIR}/lsan.supp"
# leaks related to DeePMD-kit and TensorFlow
leak:DP_NewDeepPot
leak:deepmd::AtomMap::AtomMap
EOF
fi

load "${BUILDDIR}/setup_deepmd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "deepmd"
