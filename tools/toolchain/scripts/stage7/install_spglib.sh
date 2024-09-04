#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"
spglib_ver="2.3.1"
spglib_sha256="c295dbea7d2fc9e50639aa14331fef277878c35f00ef0766e688bfbb7b17d44c"

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
    echo "==================== Installing spglib ===================="
    pkg_install_dir="${INSTALLDIR}/spglib-${spglib_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "spglib-${spglib_ver} is already installed, skipping it."
    else
      if [ -f spglib-${spglib_ver}.tar.gz ]; then
        echo "spglib-${spglib_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${spglib_sha256}" "spglib-${spglib_ver}.tar.gz"
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf spglib-${spglib_ver} "${pkg_install_dir}"
      tar -xzf spglib-${spglib_ver}.tar.gz
      cd spglib-${spglib_ver}

      mkdir build
      cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DSPGLIB_SHARED_LIBS=OFF \
        -DSPGLIB_USE_OMP=ON \
        -DSPGLIB_WITH_TESTS=OFF \
        .. > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi

    SPGLIB_CFLAGS="-I${pkg_install_dir}/include"
    SPGLIB_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    if [ -d "${pkg_install_dir}/lib64" ]; then
      ln -sf lib64 ${pkg_install_dir}/lib
      cd ${pkg_install_dir}
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding spglib from system paths ===================="
    check_command pkg-config --modversion spglib
    add_include_from_paths SPGLIB_CFLAGS "spglib.h" $INCLUDE_PATHS
    add_lib_from_paths SPGLIB_LDFLAGS "libspglib.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking spglib to user paths ===================="
    pkg_install_dir="$with_spglib"
    check_dir "$pkg_install_dir/lib"
    check_dir "$pkg_install_dir/include"
    SPGLIB_CFLAGS="-I'${pkg_install_dir}/include'"
    SPGLIB_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_spglib" != "__DONTUSE__" ]; then
  SPGLIB_LIBS="-lsymspg"
  if [ "$with_spglib" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_spglib"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_spglib"
export SPGLIB_CFLAGS="-I$pkg_install_dir/include ${SPGLIB_CFLAGS}"
export SPGLIB_LDFLAGS="${SPGLIB_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__SPGLIB"
export CP_CFLAGS="\${CP_CFLAGS} ${SPGLIB_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${SPGLIB_LDFLAGS}"
export CP_LIBS="${SPGLIB_LIBS} \${CP_LIBS}"
export LIBSPG_ROOT="$pkg_install_dir"
export LIBSPG_INCLUDE_DIR="$pkg_install_dir/include"
EOF
  cat "${BUILDDIR}/setup_spglib" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_spglib"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "spglib"
