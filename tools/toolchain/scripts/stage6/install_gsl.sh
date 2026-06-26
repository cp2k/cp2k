#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

gsl_ver="2.8"
gsl_sha256="6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190"
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
    echo "==================== Installing GSL ===================="
    pkg_install_dir="${INSTALLDIR}/gsl-${gsl_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "gsl-${gsl_ver} is already installed, skipping it."
    else
      retrieve_package "${gsl_sha256}" "gsl-${gsl_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d gsl-${gsl_ver} ] && rm -rf gsl-${gsl_ver}
      tar -xzf gsl-${gsl_ver}.tar.gz
      cd gsl-${gsl_ver}
      ./configure --prefix="${pkg_install_dir}" \
        --libdir="${pkg_install_dir}/lib" \
        --disable-shared \
        --enable-static \
        > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail_excerpt install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    GSL_CFLAGS="-I${pkg_install_dir}/include"
    GSL_LDFLAGS="-L${pkg_install_dir}/lib -Wl,-rpath,${pkg_install_dir}/lib"
    ;;
  __SYSTEM__)
    echo "==================== Finding GSL from system paths ===================="
    check_pkgconfig gsl
    pkg_install_dir="$(pkg-config --variable=prefix gsl)"
    GSL_INCLUDEDIR="$(pkg-config --variable=includedir gsl)"
    GSL_LIBDIR="$(pkg-config --variable=libdir gsl)"
    GSL_CFLAGS="-I${GSL_INCLUDEDIR}"
    GSL_LDFLAGS="-L${GSL_LIBDIR} -Wl,-rpath,${GSL_LIBDIR}"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking GSL to user paths ===================="
    pkg_install_dir="$with_gsl"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    GSL_CFLAGS="-I${pkg_install_dir}/include"
    GSL_LDFLAGS="-L${pkg_install_dir}/lib -Wl,-rpath,${pkg_install_dir}/lib"
    ;;
esac
if [ "$with_gsl" != "__DONTUSE__" ]; then
  if [ "$with_gsl" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_gsl"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_gsl"
export GSL_ROOT="${pkg_install_dir}"
export GSL_VER="${gsl_ver}"
export GSL_CFLAGS="${GSL_CFLAGS}"
export GSL_LDFLAGS="${GSL_LDFLAGS}"
EOF
  filter_setup "${BUILDDIR}/setup_gsl" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_gsl"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "gsl"
