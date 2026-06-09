#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxstream_ver="16bb679"
libxstream_sha256="6fb9f7f39fd7ae7a810bd9a4643881d8feb60ec6d62079a04ce0396f39d06df4"
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
    install_lock_file="$pkg_install_dir/install_successful"

    # libxstream depends on libxs
    if [ "$with_libxs" = "__DONTUSE__" ]; then
      report_error $LINENO "libxstream requires libxs, but libxs is disabled."
      exit 1
    fi

    if verify_checksums "${install_lock_file}"; then
      echo "libxstream-${libxstream_ver} is already installed, skipping it."
    else
      if [ -f libxstream-${libxstream_ver}.tar.gz ]; then
        echo "libxstream-${libxstream_ver}.tar.gz is found"
      else
        if ! download_pkg_from_cp2k_org "${libxstream_sha256}" "libxstream-${libxstream_ver}.tar.gz" 2> /dev/null; then
          download_pkg_from_urlpath "${libxstream_sha256}" "${libxstream_ver}" \
            https://codeload.github.com/hfp/libxstream/tar.gz \
            "libxstream-${libxstream_ver}.tar.gz"
        fi
      fi
      [ -d libxstream-${libxstream_ver} ] && rm -rf libxstream-${libxstream_ver}
      tar -xzf libxstream-${libxstream_ver}.tar.gz

      echo "Installing from scratch into ${pkg_install_dir}"
      cd libxstream-${libxstream_ver}
      # Determine LIBXS install location.
      libxs_prefix="${LIBXSROOT:-}"
      if [ -z "${libxs_prefix}" ] && command -v pkg-config > /dev/null 2>&1 && pkg-config --exists libxs; then
        libxs_prefix="$(pkg-config --variable=prefix libxs)"
      fi
      if [ -z "${libxs_prefix}" ] && [ -d "${INSTALLDIR}/libxs-${LIBXS_VER:-${libxs_ver:-unknown}}" ]; then
        libxs_prefix="${INSTALLDIR}/libxs-${LIBXS_VER:-${libxs_ver:-unknown}}"
      fi
      if [ -z "${LIBXSROOT}" ] && [ -n "${libxs_prefix}" ]; then
        LIBXSROOT="${libxs_prefix}"
      fi
      if [ -z "${LIBXSROOT}" ]; then
        report_error $LINENO "LIBXSTREAM requires LIBXSROOT. Install LIBXS first or set LIBXSROOT."
        exit 1
      fi
      mkdir build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"

      # ---- macOS: pkg-config files must not use GNU ld's "-l:libfoo.a" syntax ----
      if [[ "$(uname -s)" == "Darwin" ]]; then
        perl -pi.bak -e 's/-l:lib([A-Za-z0-9_]+)\.a\b/-l$1/g' \
          $(find "${pkg_install_dir}/lib" -name '*.pc' -type f)
      fi

    fi
    LIBXSTREAM_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXSTREAM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXStream from system paths ===================="
    check_lib -lxstream "libxstream"
    add_include_from_paths LIBXSTREAM_CFLAGS "libxstream.h" $INCLUDE_PATHS
    add_lib_from_paths LIBXSTREAM_LDFLAGS "libxstream.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBXStream to user paths ===================="
    pkg_install_dir="$with_libxstream"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    LIBXSTREAM_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXSTREAM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libxstream" != "__DONTUSE__" ]; then
  LIBXSTREAM_LIBS="-lxstream"
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
