#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

gmp_ver="6.3.0"
gmp_sha256="e56fd59d76810932a0555aa15a14b61c16bed66110d3c75cc2ac49ddaa9ab24c"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_gmp" ] && rm "${BUILDDIR}/setup_gmp"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_gmp" in
  __INSTALL__)
    echo "==================== Installing GMP ===================="
    pkg_install_dir="${INSTALLDIR}/gmp-${gmp_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "gmp-${gmp_ver} is already installed, skipping it."
    else
      retrieve_package "${gmp_sha256}" "gmp-${gmp_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d gmp-${gmp_ver} ] && rm -rf gmp-${gmp_ver}
      tar -xzf gmp-${gmp_ver}.tar.gz
      cd gmp-${gmp_ver}
      mkdir build
      cd build
      # autotools setup, out-of-source build
      ../configure CFLAGS="${CFLAGS} -std=c17" \
        --prefix="${pkg_install_dir}" \
        --libdir="${pkg_install_dir}/lib" \
        --enable-cxx=yes \
        --includedir="${pkg_install_dir}/include" \
        > configure.log 2>&1 || tail_excerpt configure.log
      make -j "$(get_nprocs)" > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage2/$(basename "${SCRIPT_NAME}")"
      cd ..
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding GMP from system paths ===================="
    check_lib -lgmp "gmp"
    check_lib -lgmpxx "gmpxx"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libgmp.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking GMP to user paths ===================="
    pkg_install_dir="${with_gmp}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_gmp" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_gmp"
export GMP_VER="${gmp_ver}"
export GMP_ROOT="${pkg_install_dir}"
EOF
  if [ "$with_gmp" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_gmp"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_gmp" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_gmp"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "gmp"
