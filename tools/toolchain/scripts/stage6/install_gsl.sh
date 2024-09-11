#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

gsl_ver="2.8"
gls_sha256="6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_gsl" ] && rm "${BUILDDIR}/setup_gsl"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_gsl" in
  __INSTALL__)
    echo "==================== Installing gsl ===================="
    pkg_install_dir="${INSTALLDIR}/gsl-${gsl_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "gsl-${gsl_ver} is already installed, skipping it."
    else
      if [ -f gsl-${gsl_ver}.tar.gz ]; then
        echo "gsl-${gsl_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${gls_sha256}" "gsl-${gsl_ver}.tar.gz"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d gsl-${gsl_ver} ] && rm -rf gsl-${gsl_ver}
      tar -xzf gsl-${gsl_ver}.tar.gz
      cd gsl-${gsl_ver}
      ./configure --prefix="${pkg_install_dir}" \
        --libdir="${pkg_install_dir}/lib" \
        --disable-shared \
        --enable-static \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi

    GSL_CFLAGS="-I'${pkg_install_dir}/include'"
    GSL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding gsl from system paths ===================="
    check_command pkg-config --modversion gsl
    add_include_from_paths GSL_CFLAGS "gsl.h" $INCLUDE_PATHS
    add_lib_from_paths GSL_LDFLAGS "libgsl.*" $LIB_PATHS
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking gsl to user paths ===================="
    pkg_install_dir="$with_gsl"
    check_dir "$pkg_install_dir/lib"
    check_dir "$pkg_install_dir/include"
    GSL_CFLAGS="-I'${pkg_install_dir}/include'"
    GSL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_gsl" != "__DONTUSE__" ]; then
  GSL_LIBS="-lgsl"
  if [ "$with_gsl" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_gsl"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
export GSL_INCLUDE_DIR="$pkg_install_dir/include"
export GSL_LIBRARY="-lgsl"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_gsl"
export GSL_VER="${gsl_ver}"
export GSL_CFLAGS="${GSL_CFLAGS}"
export GSL_LDFLAGS="${GSL_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__GSL|)"
export CP_CFLAGS="\${CP_CFLAGS} ${GSL_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${GSL_LDFLAGS}"
export GSL_LIBRARY="-lgsl"
export GSL_ROOT="${pkg_install_dir}"
export GSL_INCLUDE_DIR="$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib64/pkgconfig"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
export CP_LIBS="IF_MPI(${GSL_LIBS}|) \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_gsl" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_gsl"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "gsl"
