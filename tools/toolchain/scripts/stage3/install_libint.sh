#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

libint_ver="2.6.0"
libint_pkg="libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}.tgz"

case "$LIBINT_LMAX" in
  4)
    libint_sha256="7c8d28bfb03920936231228b79686ba0fd87ea922c267199789bc131cf21ac08"
    ;;
  5)
    libint_sha256="1cd72206afddb232bcf2179c6229fbf6e42e4ba8440e701e6aa57ff1e871e9db"
    ;;
  6)
    libint_sha256="bea76a433cd32bde280879f73b5fc8228c78b62e3ea57ace4c6d74b65910b8af"
    ;;
  7)
    libint_sha256="3bcdcc55e1dbafe38a785d4af171df8e300bb8b7775894b57186cdf35807c334"
    ;;
  *)
    report_error "Unsupported value --libint-lmax=${LIBINT_LMAX}."
    exit 1
    ;;
esac

[ -f "${BUILDDIR}/setup_libint" ] && rm "${BUILDDIR}/setup_libint"

LIBINT_CFLAGS=""
LIBINT_LDFLAGS=""
LIBINT_LIBS=""
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
        download_pkg ${DOWNLOADER_FLAGS} ${libint_sha256} \
          https://github.com/cp2k/libint-cp2k/releases/download/v${libint_ver}/${libint_pkg}
      fi

      [ -d libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX} ] && rm -rf libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}
      tar -xzf ${libint_pkg}

      echo "Installing from scratch into ${pkg_install_dir}"
      cd libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}

      # reduce debug information to level 1 since
      # level 2 (default for -g flag) leads to very large binary size
      LIBINT_CXXFLAGS="$CXXFLAGS -g1"

      # cmake build broken with libint 2.6, uncomment for libint 2.7 and above
      #cmake . -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
      #        -DCMAKE_CXX_COMPILER="$CXX" \
      #        -DLIBINT2_INSTALL_LIBDIR="${pkg_install_dir}/lib" \
      #        -DENABLE_FORTRAN=ON \
      #        -DCXXFLAGS="$LIBINT_CXXFLAGS" > configure.log 2>&1
      #cmake --build . > cmake.log 2>&1
      #cmake --build . --target install > install.log 2>&1
      ./configure --prefix=${pkg_install_dir} \
        --with-cxx="$CXX $LIBINT_CXXFLAGS" \
        --with-cxx-optflags="$LIBINT_CXXFLAGS" \
        --enable-fortran \
        --libdir="${pkg_install_dir}/lib" \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log

      if [ "${MPI_MODE}" = "intelmpi" ]; then
        # Fix bug in makefile for Fortran module
        sed -i "s/\$(CXX) \$(CXXFLAGS)/\$(FC) \$(FCFLAGS)/g" fortran/Makefile
      fi

      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log

      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi

    LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib'"
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
    LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
    LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libint" != "__DONTUSE__" ]; then
  LIBINT_LIBS="-lint2"
  if [ "$with_libint" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libint"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_libint" >> $SETUPFILE
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
fi

load "${BUILDDIR}/setup_libint"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libint"
