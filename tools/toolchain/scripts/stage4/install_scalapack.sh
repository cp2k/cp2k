#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

scalapack_ver="2.2.1"
scalapack_sha256="4aede775fdb28fa44b331875730bcd5bab130caaec225fadeccf424c8fcb55aa"
scalapack_pkg="scalapack-${scalapack_ver}.tgz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_scalapack" ] && rm "${BUILDDIR}/setup_scalapack"

SCALAPACK_CFLAGS=''
SCALAPACK_LDFLAGS=''
SCALAPACK_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_scalapack" in
  __INSTALL__)
    echo "==================== Installing ScaLAPACK ===================="
    pkg_install_dir="${INSTALLDIR}/scalapack-${scalapack_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "scalapack-${scalapack_ver} is already installed, skipping it."
    else
      require_env MATH_LIBS
      if [ -f ${scalapack_pkg} ]; then
        echo "${scalapack_pkg} is found"
      else
        download_pkg_from_cp2k_org "${scalapack_sha256}" "${scalapack_pkg}"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d scalapack-${scalapack_ver} ] && rm -rf scalapack-${scalapack_ver}
      tar -xzf ${scalapack_pkg}

      mkdir -p "scalapack-${scalapack_ver}/build"
      pushd "scalapack-${scalapack_ver}/build" > /dev/null

      cflags=""
      fflags=""
      if ("${FC}" --version | grep -q 'GNU'); then
        cflags="-fpermissive"
        fflags=$(allowed_gfortran_flags "-fallow-argument-mismatch")
      fi
      CFLAGS=${cflags} FFLAGS=${fflags} cmake -DCMAKE_FIND_ROOT_PATH="$ROOTDIR" \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_SHARED_LIBS=NO \
        -DCMAKE_BUILD_TYPE=Release .. \
        -DBUILD_TESTING=NO \
        -DSCALAPACK_BUILD_TESTS=NO \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install >> make.log 2>&1 || tail -n ${LOG_LINES} make.log

      popd > /dev/null
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    SCALAPACK_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding ScaLAPACK from system paths ===================="
    check_lib -lscalapack "ScaLAPACK"
    add_lib_from_paths SCALAPACK_LDFLAGS "libscalapack.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking ScaLAPACK to user paths ===================="
    pkg_install_dir="$with_scalapack"
    check_dir "${pkg_install_dir}/lib"
    SCALAPACK_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_scalapack" != "__DONTUSE__" ]; then
  SCALAPACK_LIBS="-lscalapack"
  cat << EOF > "${BUILDDIR}/setup_scalapack"
export SCALAPACK_VER="${scalapack_ver}"
EOF
  if [ "$with_scalapack" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_scalapack"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
export SCALAPACK_ROOT="${pkg_install_dir}"
EOF
    cat "${BUILDDIR}/setup_scalapack" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_scalapack"
export SCALAPACK_LDFLAGS="${SCALAPACK_LDFLAGS}"
export SCALAPACK_LIBS="${SCALAPACK_LIBS}"
export SCALAPACK_ROOT="${pkg_install_dir}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${SCALAPACK_LDFLAGS}|)"
export CP_LIBS="IF_MPI(-lscalapack|) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"

# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat << EOF >> ${INSTALLDIR}/lsan.supp
# leaks related to SCALAPACK
leak:pdpotrf_
EOF

load "${BUILDDIR}/setup_scalapack"
write_toolchain_env "${INSTALLDIR}"

report_timing "scalapack"
