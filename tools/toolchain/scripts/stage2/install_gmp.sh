#!/bin/bash -e

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

gmp_ver="6.3.0"
gmp_sha256="e56fd59d76810932a0555aa15a14b61c16bed66110d3c75cc2ac49ddaa9ab24c"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}"/common_vars.sh
# shellcheck disable=SC1091
source "${SCRIPT_DIR}"/tool_kit.sh
# shellcheck disable=SC1091
source "${SCRIPT_DIR}"/signal_trap.sh
# shellcheck disable=SC1091
source "${INSTALLDIR}"/toolchain.conf
# shellcheck disable=SC1091
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_gmp" ] && rm "${BUILDDIR}/setup_gmp"

GMP_CFLAGS=""
GMP_LDFLAGS=""
GMP_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
with_gmp=${with_gmp:__DONTUSE__}
case "$with_gmp" in
  __INSTALL__)
    echo "==================== Installing GMP ===================="
    pkg_install_dir="${INSTALLDIR}/gmp-${gmp_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "gmp-${gmp_ver} is already installed, skipping it."
    else
      if [ -f gmp-${gmp_ver}.tar.gz ]; then
        echo "gmp-${gmp_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${gmp_sha256}" "gmp-${gmp_ver}.tar.gz"
      fi
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
        > configure.log 2>&1 || tail -n "${LOG_LINES}" configure.log
      make -j "$(get_nprocs)" > make.log 2>&1 || tail -n "${LOG_LINES}" make.log
      make install > install.log 2>&1 || tail -n "${LOG_LINES}" install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage2/$(basename "${SCRIPT_NAME}")"
      cd ..
    fi
    GMP_CFLAGS="-I'${pkg_install_dir}/include'"
    GMP_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding GMP from system paths ===================="
    check_lib -lgmp "gmp"
    check_lib -lgmpxx "gmpxx"
    add_include_from_paths GMP_CFLAGS "gmp.h" "$INCLUDE_PATHS"
    add_lib_from_paths GSL_LDFLAGS "libgmp.*" "$LIB_PATHS"
    add_lib_from_paths GSL_LDFLAGS "libgmpxx.*" "$LIB_PATHS"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking GMP to user paths ===================="
    pkg_install_dir="$with_gmp"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    GMP_CFLAGS="-I'${pkg_install_dir}/include'"
    GMP_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_gmp" != "__DONTUSE__" ]; then
  GMP_LIBS=" -lgmp -lgmpxx"
  cat << EOF > "${BUILDDIR}/setup_gmp"
export GMP_VER="${gmp_ver}"
EOF
  if [ "$with_gmp" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_gmp"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  # TODO : Is it necessary to have a separate feature for GMP, i.e -D__GMP during compilation?
  cat << EOF >> "${BUILDDIR}/setup_gmp"
export GMP_CFLAGS="${GMP_CFLAGS}"
export GMP_LDFLAGS="${GMP_LDFLAGS}"
export GMP_LIBS="${GMP_LIBS}"
export CP_CFLAGS="\${CP_CFLAGS} ${GMP_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${GMP_LDFLAGS}"
export CP_LIBS="${GMP_LIBS} \${CP_LIBS}"
export GMP_ROOT="${pkg_install_dir}"
EOF
  cat "${BUILDDIR}/setup_gmp" >> "$SETUPFILE"
fi

load "${BUILDDIR}/setup_gmp"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "gmp"
