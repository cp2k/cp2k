#!/bin/bash -e

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libsmeagol_ver="1.2"
libsmeagol_sha256="0b76198a48c47256cd8b048f32c7752a1de4e533ed1c1f4ea097255dd712917b"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libsmeagol" ] && rm "${BUILDDIR}/setup_libsmeagol"

LIBSMEAGOL_CFLAGS=""
LIBSMEAGOL_LDFLAGS=""
LIBSMEAGOL_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libsmeagol}" in
  __INSTALL__)
    echo "==================== Installing libsmeagol ===================="
    pkg_install_dir="${INSTALLDIR}/libsmeagol-${libsmeagol_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libsmeagol-${libsmeagol_ver} is already installed, skipping it."
    else
      if [ -f libsmeagol-${libsmeagol_ver}.tar.gz ]; then
        echo "libsmeagol-${libsmeagol_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${libsmeagol_sha256}" "libsmeagol-${libsmeagol_ver}.tar.gz"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libsmeagol-${libsmeagol_ver} ] && rm -rf libsmeagol-${libsmeagol_ver}
      tar -xzf libsmeagol-${libsmeagol_ver}.tar.gz
      cd libsmeagol-${libsmeagol_ver}
      if ("${FC}" --version | grep -q 'GNU'); then
        FCFLAGS_FIX="-ffixed-form"
        FCFLAGS_FREE="-ffree-form -ffree-line-length-none"
      else
        FCFLAGS_FIX=""
        FCFLAGS_FREE=""
      fi
      make -j $(get_nprocs) \
        FC=${MPIFC} \
        FCFLAGS="-DMPI -fallow-argument-mismatch ${FCFLAGS}" \
        FCFLAGS_FIXEDFORM="${FCFLAGS_FIX}" \
        FCFLAGS_FREEFORM="${FCFLAGS_FREE}" \
        > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      # The libsmeagol makefile does not provide an install target
      [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
      mkdir ${pkg_install_dir} && cp -a lib ${pkg_install_dir}
      mkdir ${pkg_install_dir}/include && cp obj/*.mod ${pkg_install_dir}/include
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi
    LIBSMEAGOL_CFLAGS="-I${pkg_install_dir}/include"
    LIBSMEAGOL_LDFLAGS="-L${pkg_install_dir}/lib -Wl,-rpath,${pkg_install_dir}/lib"
    ;;
  __SYSTEM__)
    echo "==================== Finding libsmeagol from system paths ===================="
    echo "Please rerun the toolchain using --with-libsmeagol=<path> option"
    exit 1
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking libsmeagol to user paths ===================="
    pkg_install_dir="${with_libsmeagol}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/obj"
    LIBSMEAGOL_CFLAGS="-I'${pkg_install_dir}/obj'"
    LIBSMEAGOL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac

if [ "${with_libsmeagol}" != "__DONTUSE__" ]; then
  LIBSMEAGOL_LIBS="-l:libsmeagol.a"
  cat << EOF > "${BUILDDIR}/setup_libsmeagol"
export LIBSMEAGOL_VER="${libsmeagol_ver}"
EOF
  if [ "${with_libsmeagol}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libsmeagol"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  # libsmeagol depends on an MPI library. In its current state libsmeagol cannot be linked with a non-MPI code.
  cat << EOF >> "${BUILDDIR}/setup_libsmeagol"
export LIBSMEAGOL_CFLAGS="${LIBSMEAGOL_CFLAGS}"
export LIBSMEAGOL_LDFLAGS="${LIBSMEAGOL_LDFLAGS}"
export LIBSMEAGOL_LIBS="${LIBSMEAGOL_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__SMEAGOL|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${LIBSMEAGOL_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${LIBSMEAGOL_LDFLAGS}|)"
export CP_LIBS="IF_MPI(${LIBSMEAGOL_LIBS}|) \${CP_LIBS}"
export LIBSMEAGOL_ROOT="${pkg_install_dir}"
EOF
  cat "${BUILDDIR}/setup_libsmeagol" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_libsmeagol"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libsmeagol"
