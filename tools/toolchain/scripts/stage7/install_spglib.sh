#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"
spglib_ver="2.7.0"
spglib_sha256="b22fc9abae9716c574fbc6d55cfc53ed654a714fccc5657a26ff5d18114bd8bd"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_spglib" ] && rm "${BUILDDIR}/setup_spglib"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_spglib" in
  __INSTALL__)
    echo "==================== Installing Spglib ===================="
    pkg_install_dir="${INSTALLDIR}/spglib-${spglib_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "spglib-${spglib_ver} is already installed, skipping it."
    else
      retrieve_package "${spglib_sha256}" "spglib-${spglib_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf spglib-${spglib_ver} "${pkg_install_dir}"
      tar -xzf spglib-${spglib_ver}.tar.gz
      cd spglib-${spglib_ver}

      mkdir build
      cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DSPGLIB_SHARED_LIBS=OFF \
        -DSPGLIB_USE_OMP=ON \
        -DSPGLIB_WITH_Fortran=ON \
        -DSPGLIB_WITH_TESTS=OFF \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Spglib from system paths ===================="
    check_pkgconfig spglib
    pkg_install_dir="$(pkg-config --variable=prefix spglib)"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Spglib to user paths ===================="
    pkg_install_dir="${with_spglib}"
    SPGLIB_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && SPGLIB_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${SPGLIB_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_spglib" != "__DONTUSE__" ]; then
  if [ "$with_spglib" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_spglib"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_spglib"
export SPGLIB_ROOT="${pkg_install_dir}"
export SPGLIB_VER="${spglib_ver}"
EOF
  filter_setup "${BUILDDIR}/setup_spglib" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_spglib"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "spglib"
