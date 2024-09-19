#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libxsmm_ver="e0c4a2389afba36c453233ad7de07bd92c715bec"
libxsmm_sha256="7140650d7ce58be08af2af3f49d27641f35a131e0f712c033e5c787ac9916c9d"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libxsmm" ] && rm "${BUILDDIR}/setup_libxsmm"

LIBXSMM_CFLAGS=''
LIBXSMM_LDFLAGS=''
LIBXSMM_LIBS=''
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
        if ! download_pkg_from_cp2k_org "${libxsmm_sha256}" "libxsmm-${libxsmm_ver}.tar.gz"; then
          download_pkg_from_urlpath "${libxsmm_sha256}" "${libxsmm_ver}.tar.gz" \
            https://github.com/libxsmm/libxsmm/archive \
            "libxsmm-${libxsmm_ver}.tar.gz"
        fi
      fi
      [ -d libxsmm-${libxsmm_ver} ] && rm -rf libxsmm-${libxsmm_ver}
      tar -xzf libxsmm-${libxsmm_ver}.tar.gz

      echo "Installing from scratch into ${pkg_install_dir}"
      # note that we do not have to set -L flags to ld for the
      # linked math libraries for the libxsmm build, as for a
      # library this is not required, you just have to provide
      # the appropriate -L flags (LDFLAGS) during the linking
      # stage of building an executable that uses the libxsmm
      # library
      cd libxsmm-${libxsmm_ver}
      # Avoid an unintended (incompatible) setting of FORTRAN
      unset FORTRAN
      make -j $(get_nprocs) \
        CXX=$CXX \
        CC=$CC \
        FC=$FC \
        WRAP=0 \
        PREFIX=${pkg_install_dir} \
        > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) \
        CXX=$CXX \
        CC=$CC \
        FC=$FC \
        WRAP=0 \
        PREFIX=${pkg_install_dir} \
        install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
      mkdir ${pkg_install_dir}/lib/pkgconfig
      cp ${pkg_install_dir}/lib/*.pc ${pkg_install_dir}/lib/pkgconfig
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
EOF
  if [ "$with_libxsmm" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libxsmm"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
EOF
    cat "${BUILDDIR}/setup_libxsmm" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_libxsmm"
export LIBXSMM_CFLAGS="${LIBXSMM_CFLAGS}"
export LIBXSMM_LDFLAGS="${LIBXSMM_LDFLAGS}"
export LIBXSMM_LIBS="${LIBXSMM_LIBS}"
export CP_DFLAGS="-D__LIBXSMM \${CP_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBXSMM_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBXSMM_LDFLAGS}"
export CP_LIBS="\${LIBXSMM_LIBS} \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"

load "${BUILDDIR}/setup_libxsmm"
write_toolchain_env "${INSTALLDIR}"

report_timing "libxsmm"
