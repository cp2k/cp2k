#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libgrpp_ver="20231215"
libgrpp_sha="3b0f55795a0c2699b81f403b1c81c56e26332f780a87655410745a0ccb51ef2f"
libgrpp_pkg="libgrpp-main-${libgrpp_ver}.zip"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libgrpp" ] && rm "${BUILDDIR}/setup_libgrpp"

LIBGRPP_CFLAGS=""
LIBGRPP_LDFLAGS=""
LIBGRPP_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libgrpp}" in
  __INSTALL__)
    echo "==================== Installing LIBGRPP ===================="
    pkg_install_dir="${INSTALLDIR}/libgrpp-main-${libgrpp_ver}"
    install_lock_file="$pkg_install_dir/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "libgrpp-main-${libgrpp_ver} is already installed, skipping it."
    else
      if [ -f ${libgrpp_pkg} ]; then
        echo "${libgrpp_pkg} is found"
      else
        download_pkg_from_cp2k_org "${libgrpp_sha}" "${libgrpp_pkg}"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libgrpp-main ] && rm -rf libgrpp-main
      unzip -qq ${libgrpp_pkg}
      cd libgrpp-main
      mkdir build
      cd build
      CC=${CC} FC=${FC} cmake .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      cd ..

      install -d "${pkg_install_dir}/lib" >> install.log 2>&1
      install -d "${pkg_install_dir}/include" >> install.log 2>&1
      install -m 644 build/libgrpp/liblibgrpp.a "${pkg_install_dir}/lib" >> install.log 2>&1
      install -m 644 build/libgrpp.mod "${pkg_install_dir}/include" >> install.log 2>&1

      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi

    LIBGRPP_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBGRPP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding libgrpp from system paths ===================="
    check_lib -llibgrpp "libgrpp"
    add_include_from_paths -p LIBGRPP_CFLAGS "libgrpp" $INCLUDE_PATHS
    add_lib_from_paths LIBGRPP_LDFLAGS "liblibgrpp.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking libgrpp to user paths ===================="
    pkg_install_dir="$with_libgrpp"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    LIBGRPP_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBGRPP_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libgrpp" != "__DONTUSE__" ]; then
  LIBGRPP_LIBS="-llibgrpp"
  if [ "$with_libgrpp" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libgrpp"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export LIBGRPP_ROOT="${pkg_install_dir}"
EOF
    cat "${BUILDDIR}/setup_libgrpp" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_libgrpp"
export LIBGRPP_CFLAGS="${LIBGRPP_CFLAGS}"
export LIBGRPP_LDFLAGS="${LIBGRPP_LDFLAGS}"
export LIBGRPP_LIBS="${LIBGRPP_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBGRPP"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBGRPP_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBGRPP_LDFLAGS}"
export CP_LIBS="${LIBGRPP_LIBS} \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_libgrpp"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libgrpp"
