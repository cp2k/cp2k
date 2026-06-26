#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libvdwxc_ver="0.5.0"
libvdwxc_sha256="29fb70efd58aff51524d2172a87e8f88e760b696b0ddb9aa5878432bdffa3c2f"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libvdwxc" ] && rm "${BUILDDIR}/setup_libvdwxc"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libvdwxc" in
  __INSTALL__)
    require_env FFTW3_INCLUDES
    require_env FFTW3_LIBS

    echo "==================== Installing libvdwxc ===================="
    pkg_install_dir="${INSTALLDIR}/libvdwxc-${libvdwxc_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libvdwxc-${libvdwxc_ver} is already installed, skipping it."
    else
      retrieve_package "${libvdwxc_sha256}" "libvdwxc-${libvdwxc_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libvdwxc-${libvdwxc_ver} ] && rm -rf libvdwxc-${libvdwxc_ver}
      tar -xzf libvdwxc-${libvdwxc_ver}.tar.gz
      cd libvdwxc-${libvdwxc_ver}

      if [ "$(uname -s)" = "Darwin" ]; then
        LDFLAGS="${LDFLAGS} -ld_classic"
      fi

      if [ "${MPI_MODE}" = "no" ]; then
        # compile libvdwxc without mpi support since fftw (or mkl) do not have mpi support activated
        ./configure \
          CC="${CC}" CFLAGS="${CFLAGS} -fpermissive" \
          FC="${FC}" \
          FFTW3_INCLUDES="${FFTW3_INCLUDES}" \
          FFTW3_LIBS="$(resolve_string "${FFTW3_LIBS}" "MPI")" \
          --prefix="${pkg_install_dir}" \
          --libdir="${pkg_install_dir}/lib" \
          --disable-shared \
          > configure.log 2>&1 || tail_excerpt configure.log
      else
        ./configure \
          CC="${MPICC}" CFLAGS="${CFLAGS} -fpermissive" \
          FC="${MPIFC}" \
          FFTW3_INCLUDES="${FFTW3_INCLUDES}" \
          FFTW3_LIBS="$(resolve_string "${FFTW3_LIBS}" "MPI")" \
          --prefix="${pkg_install_dir}" \
          --libdir="${pkg_install_dir}/lib" \
          --disable-shared \
          > configure.log 2>&1 || tail_excerpt configure.log
      fi
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding libvdwxc from system paths ===================="
    check_lib -lvdwxc "libvdwxc"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libvdwxc.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking libvdwxc to user paths ===================="
    pkg_install_dir="${with_libvdwxc}"
    VDWXC_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && VDWXC_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${VDWXC_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_libvdwxc" != "__DONTUSE__" ]; then
  if [ "$with_libvdwxc" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libvdwxc"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libvdwxc"
export LIBVDWXC_VER="${libvdwxc_ver}"
export LIBVDWXC_ROOT="${pkg_install_dir}"
export VDWXC_ROOT="${pkg_install_dir}"
EOF
  filter_setup "${BUILDDIR}/setup_libvdwxc" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libvdwxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libvdwxc"
