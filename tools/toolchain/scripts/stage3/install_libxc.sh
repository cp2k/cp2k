#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxc_ver="7.1.2"
libxc_sha256="3915fac94416e4c415534223ea492ad2663f928acf27e98662c861b094a6c306"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxc" ] && rm "${BUILDDIR}/setup_libxc"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libxc" in
  __INSTALL__)
    echo "==================== Installing LIBXC ===================="
    pkg_install_dir="${INSTALLDIR}/libxc-${libxc_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libxc-${libxc_ver} is already installed, skipping it."
    else
      retrieve_package "${libxc_sha256}" "libxc-${libxc_ver}.tar.bz2"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libxc-${libxc_ver} ] && rm -rf libxc-${libxc_ver}
      tar -xjf libxc-${libxc_ver}.tar.bz2
      cd libxc-${libxc_ver}
      mkdir build
      cd build
      if [ "${with_gcc}" != "__DONTUSE__" ] &&
        [ "${with_intel}" = "__DONTUSE__" ] && [ "${with_amd}" = "__DONTUSE__" ]; then
        # Turn off variable tracking
        LIBXC_CFLAGS="${CFLAGS} -fno-var-tracking"
      else
        LIBXC_CFLAGS=""
      fi
      # CP2K make use of third derivatives in libxc
      CFLAGS="${LIBXC_CFLAGS}" cmake \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DBUILD_TESTING=OFF \
        -DENABLE_FORTRAN=ON \
        -DMAXORDER=3 \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
      cd ..
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBXC from system paths ===================="
    check_lib -lxcf03 "libxc"
    check_lib -lxc "libxc"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libxc*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBXC to user paths ===================="
    pkg_install_dir="${with_libxc}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_libxc" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libxc"
export LIBXC_ROOT="${pkg_install_dir}"
export LIBXC_VER="${libxc_ver}"
EOF
  if [ "$with_libxc" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxc"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libxc" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libxc"
