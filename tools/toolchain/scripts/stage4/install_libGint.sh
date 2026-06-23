#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"
libgint_ver="v1"
libgint_sha256=cc0dfeb6022ebfe0c3028a131045ff49b4c1005ad9cb77c5736fb3dac045b192

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libGint" ] && rm "${BUILDDIR}/setup_libGint"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libgint" in
  __INSTALL__)
    echo "==================== Installing libGint ===================="
    pkg_install_dir="${INSTALLDIR}/libGint-${libgint_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libGint-${libgint_ver} is already installed, skipping it."
    else
      retrieve_package "${libgint_sha256}" "libGint-${libgint_ver}.tar.gz"
      # The tar called libGint-v1 expands into libGint-release_v1
      [ -d libGint-release_${libgint_ver} ] && rm -rf libGint-release_${libgint_ver}
      tar -xzf libGint-${libgint_ver}.tar.gz
      echo "Installing from scratch into ${pkg_install_dir}"
      cd libGint-release_${libgint_ver}

      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install PREFIX="${pkg_install_dir}" > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding libGint from system paths ===================="
    check_lib -lGint "libGint"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libxs.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking libGint to user paths ===================="
    echo "$with_libgint"
    pkg_install_dir="$with_libgint"
    check_dir "$pkg_install_dir/lib"
    check_dir "$pkg_install_dir/include"
    ;;
esac

if [ "$with_libgint" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libGint"
export LIBGINT_VER="${libgint_ver}"
EOF
  if [ "$with_libgint" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libGint"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libGint" "${SETUPFILE}"
fi
cd "${ROOTDIR}"

load "${BUILDDIR}/setup_libGint"
write_toolchain_env "${INSTALLDIR}"

report_timing "libGint"
