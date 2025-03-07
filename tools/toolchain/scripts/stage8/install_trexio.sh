#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

trexio_ver="2.5.0"
trexio_sha256="7bf7e0021467530b4946fb3f6ee39f393e2f4bd65a6f4debaec774120c29e4ee"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_trexio" ] && rm "${BUILDDIR}/setup_trexio"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_trexio" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing TREXIO ===================="
    require_env HDF5_LIBS
    require_env HDF5_CFLAGS
    require_env HDF5_LDFLAGS

    pkg_install_dir="${INSTALLDIR}/trexio-${trexio_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "trexio-${trexio_ver} is already installed, skipping it."
    else
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d trexio-${trexio_ver} ] && rm -rf trexio-${trexio_ver}

      if [ -f trexio-${trexio_ver}.tar.gz ]; then
        echo "trexio_${trexio_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${trexio_sha256}" "trexio-${trexio_ver}.tar.gz"
      fi

      tar -xzf trexio-${trexio_ver}.tar.gz
      cd trexio-${trexio_ver}

      ./configure --prefix="${pkg_install_dir}" --libdir="${pkg_install_dir}/lib" > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log

      make -j $(get_nprocs) >> make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    TREXIO_CFLAGS="-I${pkg_install_dir}/include"
    TREXIO_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding trexio from system paths ===================="
    require_env HDF5_LIBS
    require_env HDF5_CFLAGS
    require_env HDF5_LDFLAGS
    check_lib -ltrexio "trexio"
    add_include_from_paths trexio_CFLAGS "trexio*" $INCLUDE_PATHS
    add_lib_from_paths trexio_LDFLAGS "libtrexio.*" $LIB_PATHS
    ;;
  *)
    echo "==================== Linking TREXIO to user paths ===================="
    pkg_install_dir="$with_trexio"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    TREXIO_CFLAGS="-I'${pkg_install_dir}/include'"
    TREXIO_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_trexio" != "__DONTUSE__" ]; then
  TREXIO_LIBS="-l:libtrexio.a"
  cat << EOF > "${BUILDDIR}/setup_trexio"
export TREXIO_VER="${trexio_ver}"
EOF
  if [ "$with_trexio" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_trexio"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_trexio"
export TREXIO_CFLAGS="${TREXIO_CFLAGS}"
export TREXIO_LDFLAGS="${TREXIO_LDFLAGS}"
export TREXIO_LIB="${TREXIO_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__TREXIO"
export CP_CFLAGS="\${CP_CFLAGS} ${TREXIO_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${TREXIO_LDFLAGS}"
export CP_LIBS="${TREXIO_LIBS} \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_trexio" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_trexio"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "trexio"
