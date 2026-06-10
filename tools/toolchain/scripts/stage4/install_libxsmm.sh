#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxsmm_ver="cdeedf7"
libxsmm_sha256="4190ce3c3ede3721b3b23e0985a618805131befe42b8deace6ed479e4fe21004"
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
      if [ -f libxsmm-${libxsmm_ver}.tar.gz ]; then
        echo "libxsmm-${libxsmm_ver}.tar.gz is found"
      else
        if ! download_pkg_from_cp2k_org "${libxsmm_sha256}" "libxsmm-${libxsmm_ver}.tar.gz" 2> /dev/null; then
          download_pkg_from_urlpath "${libxsmm_sha256}" "${libxsmm_ver}" \
            https://codeload.github.com/libxsmm/libxsmm/tar.gz \
            "libxsmm-${libxsmm_ver}.tar.gz"
        fi
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libxsmm-${libxsmm_ver} ] && rm -rf libxsmm-${libxsmm_ver}
      tar -xzf libxsmm-${libxsmm_ver}.tar.gz
      cd libxsmm-${libxsmm_ver}
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
    LIBXSMM_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXSMM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding Libxsmm from system paths ===================="
    check_lib -lxsmm "libxsmm"
    check_lib -lxsmmf "libxsmm"
    check_lib -lxsmmext "libxsmm"
    add_include_from_paths LIBXSMM_CFLAGS "libxsmm.h" $INCLUDE_PATHS
    add_lib_from_paths LIBXSMM_LDFLAGS "libxsmm.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Libxsmm to user paths ===================="
    pkg_install_dir="$with_libxsmm"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    LIBXSMM_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXSMM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libxsmm" != "__DONTUSE__" ]; then
  LIBXSMM_LIBS="-lxsmmf -lxsmmext -lxsmm -ldl -lpthread"
  cat << EOF > "${BUILDDIR}/setup_libxsmm"
export LIBXSMM_VER="${libxsmm_ver}"
export LIBXSMM_ROOT="${pkg_install_dir:-}"
EOF
  if [ "$with_libxsmm" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxsmm"
prepend_path PATH "${pkg_install_dir}/bin"
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
