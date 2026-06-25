#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxsmm_ver="db07b74"
libxsmm_sha256="741c1ca39deee4bd01fc62c076b68a5623bb5e1aeaea0368fb65d2fdf7bdbed7"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxsmm" ] && rm "${BUILDDIR}/setup_libxsmm"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxsmm" in
  __INSTALL__)
    echo "==================== Installing Libxsmm ===================="
    if [[ ("$OPENBLAS_ARCH" != "x86_64") && ("$OPENBLAS_ARCH" != "arm64") ]]; then
      report_warning $LINENO "libxsmm is not supported on arch ${OPENBLAS_ARCH}"
      cat << EOF > "${BUILDDIR}/setup_libxsmm"
with_libxsmm="__DONTUSE__"
EOF
      exit 0
    fi
    pkg_install_dir="${INSTALLDIR}/libxsmm-${libxsmm_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libxsmm-${libxsmm_ver} is already installed, skipping it."
    else
      retrieve_package "${libxsmm_sha256}" "libxsmm-${libxsmm_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libxsmm-${libxsmm_ver} ] && rm -rf libxsmm-${libxsmm_ver}
      tar -xzf libxsmm-${libxsmm_ver}.tar.gz
      cd libxsmm-${libxsmm_ver}
      mkdir build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DLIBXSMM_FORTRAN=ON \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Libxsmm from system paths ===================="
    check_lib -lxsmm "libxsmm"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libxsmm.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Libxsmm to user paths ===================="
    pkg_install_dir="$with_libxsmm"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    ;;
esac
if [ "$with_libxsmm" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libxsmm"
export LIBXSMM_VER="${libxsmm_ver}"
export LIBXSMM_ROOT="${pkg_install_dir:-}"
EOF
  if [ "$with_libxsmm" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxsmm"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libxsmm" "${SETUPFILE}"
fi
cd "${ROOTDIR}"

load "${BUILDDIR}/setup_libxsmm"
write_toolchain_env "${INSTALLDIR}"

report_timing "libxsmm"
