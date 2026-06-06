#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxstream_ver="8015f5461e5d1fee08024ddf20b38ea4400fbc24"
libxstream_sha256="f2c41d2c0cbe1fd6d8001d78834caf1e978acdec61e4d873373684be6e98d981"
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
          download_pkg_from_urlpath "${libxstream_sha256}" "${libxstream_ver}.tar.gz" \
            https://github.com/hfp/libxstream/archive \
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
      make -j $(get_nprocs) \
        CXX=$CXX \
        CC=$CC \
        LIBXSROOT="${LIBXSROOT}" \
        PREFIX=${pkg_install_dir} \
        > make.log 2>&1 || tail_excerpt make.log
      make -j $(get_nprocs) \
        CXX=$CXX \
        CC=$CC \
        LIBXSROOT="${LIBXSROOT}" \
        PREFIX=${pkg_install_dir} \
        install > install.log 2>&1 || tail_excerpt install.log
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
export LIBXSTREAMROOT="${pkg_install_dir}"
EOF
  if [ "$with_libxstream" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxstream"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libxstream"
export LIBXSTREAM_CFLAGS="${LIBXSTREAM_CFLAGS}"
export LIBXSTREAM_LDFLAGS="${LIBXSTREAM_LDFLAGS}"
export LIBXSTREAM_LIBS="${LIBXSTREAM_LIBS}"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBXSTREAM_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBXSTREAM_LDFLAGS}"
export CP_LIBS="\${LIBXSTREAM_LIBS} \${CP_LIBS}"
EOF
  filter_setup "${BUILDDIR}/setup_libxstream" "${SETUPFILE}"
fi
cd "${ROOTDIR}"

load "${BUILDDIR}/setup_libxstream"
write_toolchain_env "${INSTALLDIR}"

report_timing "libxstream"
