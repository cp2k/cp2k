#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxs_ver="81914e7"
libxs_sha256="afa98ef4f3fab9c5500fb3a5a3179aabf1ff29013da8b201f164afd2215ff996"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxs" ] && rm "${BUILDDIR}/setup_libxs"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxs" in
  __INSTALL__)
    echo "==================== Installing LIBXS ===================="
    pkg_install_dir="${INSTALLDIR}/libxs-${libxs_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libxs-${libxs_ver} is already installed, skipping it."
    else
      if [ -f libxs-${libxs_ver}.tar.gz ]; then
        echo "libxs-${libxs_ver}.tar.gz is found"
      else
        download_pkg_from_urlpath "${libxs_sha256}" "${libxs_ver}" \
          https://codeload.github.com/hfp/libxs/tar.gz \
          "libxs-${libxs_ver}.tar.gz"
      fi
      [ -d libxs-${libxs_ver} ] && rm -rf libxs-${libxs_ver}
      tar -xzf libxs-${libxs_ver}.tar.gz

      echo "Installing from scratch into ${pkg_install_dir}"
      cd libxs-${libxs_ver}
      mkdir build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DLIBXS_FORTRAN=ON \
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
    LIBXS_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXS from system paths ===================="
    check_lib -lxs "libxs"
    add_include_from_paths LIBXS_CFLAGS "libxs.h" $INCLUDE_PATHS
    add_lib_from_paths LIBXS_LDFLAGS "libxs.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBXS to user paths ===================="
    pkg_install_dir="$with_libxs"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    LIBXS_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libxs" != "__DONTUSE__" ]; then
  LIBXS_LIBS="-lxs -ldl -lpthread"
  cat << EOF > "${BUILDDIR}/setup_libxs"
export LIBXS_VER="${libxs_ver}"
export LIBXS_ROOT="${pkg_install_dir}"
EOF
  if [ "$with_libxs" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxs"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libxs" "${SETUPFILE}"
fi
cd "${ROOTDIR}"

load "${BUILDDIR}/setup_libxs"
write_toolchain_env "${INSTALLDIR}"

report_timing "libxs"
