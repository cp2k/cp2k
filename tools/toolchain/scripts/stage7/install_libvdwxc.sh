#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"
libvdwxc_ver="0.4.0"
libvdwxc_sha256="3524feb5bb2be86b4688f71653502146b181e66f3f75b8bdaf23dd1ae4a56b33"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libvdwxc" ] && rm "${BUILDDIR}/setup_libvdwxc"

if [ "$MPI_MODE" = "no" ] && [ $with_sirius = "__FALSE__" ]; then
  report_warning $LINENO "MPI and SIRIUS are disabled, skipping libvdwxc installation"
  exit 0
fi

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libvdwxc" in
  __INSTALL__)
    require_env FFTW3_INCLUDES
    require_env FFTW3_LIBS
    require_env FFTW_LDFLAGS
    require_env FFTW_LIBS
    require_env FFTW_CFLAGS
    require_env MPI_CFLAGS
    require_env MPI_LDFLAGS
    require_env MPI_LIBS

    echo "==================== Installing libvdwxc ===================="
    pkg_install_dir="${INSTALLDIR}/libvdwxc-${libvdwxc_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libvdwxc-${libvdwxc_ver} is already installed, skipping it."
    else
      if [ -f libvdwxc-${libvdwxc_ver}.tar.gz ]; then
        echo "libvdwxc-${libvdwxc_ver}.tar.gz is found"
      else
        # do not remove this. They do not publish official version often
        download_pkg ${DOWNLOADER_FLAGS} ${libvdwxc_sha256} \
          "https://www.cp2k.org/static/downloads/libvdwxc-${libvdwxc_ver}.tar.gz"
      fi

      for patch in "${patches[@]}"; do
        fname="${patch##*/}"
        if [ -f "${fname}" ]; then
          echo "${fname} is found"
        else
          # parallel build patch
          download_pkg ${DOWNLOADER_FLAGS} "${patch}"
        fi
      done

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libvdwxc-${libvdwxc_ver} ] && rm -rf libvdwxc-${libvdwxc_ver}
      tar -xzf libvdwxc-${libvdwxc_ver}.tar.gz
      cd libvdwxc-${libvdwxc_ver}

      for patch in "${patches[@]}"; do
        patch -p1 < ../"${patch##*/}"
      done
      # unset MPICC MPICXX MPIF90 MPIFC MPIF77
      if [ "$MPI_MODE" = "no" ]; then
        # compile libvdwxc without mpi support since fftw (or mkl) do not have mpi support activated
        ./configure \
          --prefix="${pkg_install_dir}" \
          --libdir="${pkg_install_dir}/lib" \
          --disable-shared \
          --without-mpi \
          >> configure.log 2>&1
      else
        ./configure \
          FFTW3_INCLUDES="-I${FFTW_ROOT}/include" \
          FFTW3_LIBS="-L${FFTW_ROOT}/lib" \
          --prefix="${pkg_install_dir}" \
          --libdir="${pkg_install_dir}/lib" \
          > configure.log 2>&1
      fi
      make -j $(get_nprocs) > compile.log 2>&1
      make install > compile.log 2>&1
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi

    LIBVDWXC_CFLAGS="-I${pkg_install_dir}/include"
    LIBVDWXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding libvdwxc from system paths ===================="
    check_command pkg-config --modversion libvdwxc
    add_include_from_paths LIBVDWXC_CFLAGS "vdwxc.h" $INCLUDE_PATHS
    add_lib_from_paths LIBVDWXC_LDFLAGS "libvdwxc*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking libvdwxc to user paths ===================="
    pkg_install_dir="$with_libvdwxc"
    check_dir "$pkg_install_dir/lib"
    check_dir "$pkg_install_dir/include"
    LIBVDWXC_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBVDWXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libvdwxc" != "__DONTUSE__" ]; then
  LIBVDWXC_LIBS="-lvdwxc"
  if [ "$with_libvdwxc" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libvdwxc"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libvdwxc"
export LIBVDWXC_CFLAGS="-I$pkg_install_dir/include ${LIBVDWXC_CFLAGS}"
export LIBVDWXC_LDFLAGS="${LIBVDWXC_LDFLAGS}"
export LIBVDWXC_LIBS="${LIBVDWXC_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__LIBVDWXC|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${LIBVDWXC_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${LIBVDWXC_LDFLAGS}|)"
export CP_LIBS="IF_MPI(${LIBVDWXC_LIBS}|) \${CP_LIBS}"
export PKG_CONFIG_PATH="$pkg_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH"
export VDWXC_DIR="$pkg_install_dir"
EOF
  cat "${BUILDDIR}/setup_libvdwxc" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_libvdwxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libvdwxc"
