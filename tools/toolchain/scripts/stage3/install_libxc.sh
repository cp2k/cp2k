#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxc_ver="7.0.0"
libxc_sha256="e9ae69f8966d8de6b7585abd9fab588794ada1fab8f689337959a35abbf9527d"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxc" ] && rm "${BUILDDIR}/setup_libxc"

LIBXC_CFLAGS=""
LIBXC_LDFLAGS=""
LIBXC_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxc" in
  __INSTALL__)
    echo "==================== Installing LIBXC ===================="
    pkg_install_dir="${INSTALLDIR}/libxc-${libxc_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libxc-${libxc_ver} is already installed, skipping it."
    else
      if [ -f libxc-${libxc_ver}.tar.bz2 ]; then
        echo "libxc-${libxc_ver}.tar.bz2 is found"
      else
        download_pkg_from_cp2k_org "${libxc_sha256}" "libxc-${libxc_ver}.tar.bz2"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libxc-${libxc_ver} ] && rm -rf libxc-${libxc_ver}
      tar -xjf libxc-${libxc_ver}.tar.bz2
      cd libxc-${libxc_ver}
      mkdir build
      cd build
      # CP2K does not make use of fourth derivatives, so skip their compilation with -DDISABLE_KXC=OFF
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_TESTING=OFF \
        -DENABLE_FORTRAN=ON \
        -DDISABLE_KXC=OFF \
        .. > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
      cd ..
    fi
    LIBXC_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXC from system paths ===================="
    check_lib -lxcf03 "libxc"
    check_lib -lxc "libxc"
    add_include_from_paths LIBXC_CFLAGS "xc.h" $INCLUDE_PATHS
    add_lib_from_paths LIBXC_LDFLAGS "libxc.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBXC to user paths ===================="
    pkg_install_dir="$with_libxc"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    LIBXC_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libxc" != "__DONTUSE__" ]; then
  LIBXC_LIBS="-lxcf03 -lxc"
  cat << EOF > "${BUILDDIR}/setup_libxc"
export LIBXC_VER="${libxc_ver}"
EOF
  if [ "$with_libxc" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxc"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
    cat "${BUILDDIR}/setup_libxc" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_libxc"
export LIBXC_CFLAGS="${LIBXC_CFLAGS}"
export LIBXC_LDFLAGS="${LIBXC_LDFLAGS}"
export LIBXC_LIBS="${LIBXC_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBXC"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBXC_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBXC_LDFLAGS}"
export CP_LIBS="${LIBXC_LIBS} \${CP_LIBS}"
export LIBXC_ROOT="${pkg_install_dir}"
EOF
fi

load "${BUILDDIR}/setup_libxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libxc"
