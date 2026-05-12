#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libfci_ver="0.1.0"
libfci_sha256="63e8bd632e55f4b99dbaa3e905cdcd3c5dbc17aefdca1391b855f34a792bb29a"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libfci" ] && rm "${BUILDDIR}/setup_libfci"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libfci" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing libfci ===================="
    require_env MATH_LIBS

    pkg_install_dir="${INSTALLDIR}/libfci-${libfci_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libfci-${libfci_ver} is already installed, skipping it."
    else
      retrieve_package "${libfci_sha256}" "libfci-${libfci_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d "libfci-${libfci_ver}" ] && rm -rf "libfci-${libfci_ver}"
      mkdir "libfci-${libfci_ver}"
      tar -xzf "libfci-${libfci_ver}.tar.gz" -C "libfci-${libfci_ver}" --strip-components=1

      mkdir "libfci-${libfci_ver}/build"
      cd "libfci-${libfci_ver}/build"
      cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        .. > cmake.log 2>&1 || tail_excerpt cmake.log
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . > build.log 2>&1 || tail_excerpt build.log
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . --target install > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename "${SCRIPT_NAME}")"
      cd ../..
    fi
    LIBFCI_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBFCI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding libfci from system paths ===================="
    check_lib -lfci "libfci"
    add_include_from_paths LIBFCI_CFLAGS "libfci.h" $INCLUDE_PATHS
    add_lib_from_paths LIBFCI_LDFLAGS "libfci.*" $LIB_PATHS
    ;;
  *)
    echo "==================== Linking libfci to user paths ===================="
    pkg_install_dir="$with_libfci"
    LIBFCI_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBFCI_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${LIBFCI_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    LIBFCI_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBFCI_LDFLAGS="-L'${LIBFCI_LIBDIR}' -Wl,-rpath,'${LIBFCI_LIBDIR}'"
    ;;
esac
if [ "$with_libfci" != "__DONTUSE__" ]; then
  LIBFCI_LIBS="-lfci -lstdc++"
  cat << EOF > "${BUILDDIR}/setup_libfci"
export LIBFCI_VER="${libfci_ver}"
EOF
  if [ "$with_libfci" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libfci"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libfci"
export LIBFCI_CFLAGS="${LIBFCI_CFLAGS}"
export LIBFCI_LDFLAGS="${LIBFCI_LDFLAGS}"
export LIBFCI_LIBS="${LIBFCI_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBFCI"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBFCI_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBFCI_LDFLAGS}"
export CP_LIBS="${LIBFCI_LIBS} \${CP_LIBS}"
export LIBFCI_ROOT="${pkg_install_dir}"
EOF
  filter_setup "${BUILDDIR}/setup_libfci" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libfci"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libfci"
