#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

libint_ver="2.13.1"
libint_pkg="libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}.tar.xz"

case "$LIBINT_LMAX" in
  4)
    libint_sha256="97bcc62de2c33dda9e96fc52ca47a7ab17fdfa200d15417354165d5579d6932a"
    ;;
  5)
    libint_sha256="988e35e2cad5a3e28974d33185135c57fa9800a336e33052c0579f5ebb5cf354"
    ;;
  6)
    libint_sha256="eae01c15a0cf56940ba21751c09aa5be4d35eeed454cd54b629babcce243aa71"
    ;;
  7)
    libint_sha256="1f43901e9528efa06f3ab7d1b68a219d06ece2bca0d4083750f174b70fdaf64a"
    ;;
  *)
    report_error "Unsupported value --libint-lmax=${LIBINT_LMAX}."
    exit 1
    ;;
esac

[ -f "${BUILDDIR}/setup_libint" ] && rm "${BUILDDIR}/setup_libint"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libint" in
  __INSTALL__)
    echo "==================== Installing LIBINT ===================="
    pkg_install_dir="${INSTALLDIR}/libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libint-${libint_ver} is already installed, skipping it."
    else
      if [ -f ${libint_pkg} ]; then
        echo "${libint_pkg} is found"
      else
        #download_pkg_from_cp2k_org "${libint_sha256}" "${libint_pkg}"
        download_pkg_from_urlpath "${libint_sha256}" "${libint_pkg}" https://github.com/Growl1234/RTDProject/releases/download/libint-cp2k
      fi

      [ -d libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX} ] && rm -rf libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}
      tar -xJf ${libint_pkg}

      echo "Installing from scratch into ${pkg_install_dir}"
      cd libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}

      # reduce debug information to level 1 since
      # level 2 (default for -g flag) leads to very large binary size
      LIBINT_CXXFLAGS="$CXXFLAGS -g1"
      LIBINT_FCFLAGS="$FCFLAGS -g1"

      mkdir build
      cd build
      CXXFLAGS="$LIBINT_CXXFLAGS" \
        FCFLAGS="$LIBINT_FCFLAGS" cmake .. \
        -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DLIBINT2_INSTALL_LIBDIR="${pkg_install_dir}/lib" \
        -DLIBINT2_ENABLE_FORTRAN=ON \
        > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log

      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi

    LIBINT_CFLAGS="-I${pkg_install_dir}/include"
    LIBINT_LDFLAGS="-L${pkg_install_dir}/lib"
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBINT from system paths ===================="
    check_lib -lint2 "libint"
    add_include_from_paths -p LIBINT_CFLAGS "libint" $INCLUDE_PATHS
    add_lib_from_paths LIBINT_LDFLAGS "libint2.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBINT to user paths ===================="
    pkg_install_dir="$with_libint"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    LIBINT_CFLAGS="-I${pkg_install_dir}/include"
    LIBINT_LDFLAGS="-L${pkg_install_dir}/lib"
    ;;
esac
if [ "$with_libint" != "__DONTUSE__" ]; then
  LIBINT_LIBS="-lint2"
  cat << EOF > "${BUILDDIR}/setup_libint"
export LIBINT_VER="${libint_ver}"
EOF
  if [ "$with_libint" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libint"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export LIBINT2_ROOT="${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_libint"
export LIBINT_CFLAGS="${LIBINT_CFLAGS}"
export LIBINT_LDFLAGS="${LIBINT_LDFLAGS}"
export LIBINT_LIBS="${LIBINT_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBINT"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBINT_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBINT_LDFLAGS}"
export CP_LIBS="${LIBINT_LIBS} \${CP_LIBS}"
EOF
  filter_setup "${BUILDDIR}/setup_libint" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libint"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libint"
