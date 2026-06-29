#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxstream_ver="1.0.0"
libxstream_sha256="44a2823b12eb58b5eaf97649244b93dcf921597ceabc718053cd28e5f59260e3"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxstream" ] && rm "${BUILDDIR}/setup_libxstream"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxstream" in
  __INSTALL__)
    echo "==================== Installing LIBXStream ===================="
    pkg_install_dir="${INSTALLDIR}/libxstream-${libxstream_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    # libxstream depends on libxs
    if [ "$with_libxs" = "__DONTUSE__" ]; then
      report_error $LINENO "libxstream requires libxs, but libxs is disabled."
      exit 1
    fi

    if verify_checksums "${install_lock_file}"; then
      echo "libxstream-${libxstream_ver} is already installed, skipping it."
    else
      retrieve_package "${libxstream_sha256}" "libxstream-${libxstream_ver}.tar.gz"
      [ -d libxstream-${libxstream_ver} ] && rm -rf libxstream-${libxstream_ver}
      tar -xzf libxstream-${libxstream_ver}.tar.gz

      echo "Installing from scratch into ${pkg_install_dir}"
      cd libxstream-${libxstream_ver}
      mkdir build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXStream from system paths ===================="
    check_lib -lxstream "libxstream"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libxstream.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBXStream to user paths ===================="
    pkg_install_dir="$with_libxstream"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    ;;
esac
if [ "$with_libxstream" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libxstream"
export LIBXSTREAM_VER="${libxstream_ver}"
export LIBXSTREAM_ROOT="${pkg_install_dir}"
EOF
  if [ "$with_libxstream" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxstream"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libxstream" "${SETUPFILE}"
fi
cd "${ROOTDIR}"

load "${BUILDDIR}/setup_libxstream"
write_toolchain_env "${INSTALLDIR}"

report_timing "libxstream"
