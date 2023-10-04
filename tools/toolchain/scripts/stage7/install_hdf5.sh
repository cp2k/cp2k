#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

hdf5_ver="1.14.2"
hdf5_sha256="ea3c5e257ef322af5e77fc1e52ead3ad6bf3bb4ac06480dd17ee3900d7a24cfb"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_hdf5" ] && rm "${BUILDDIR}/setup_hdf5"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_hdf5" in
  __INSTALL__)
    echo "==================== Installing hdf5 ===================="
    pkg_install_dir="${INSTALLDIR}/hdf5-${hdf5_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "hdf5-${hdf5_ver} is already installed, skipping it."
    else
      if [ -f hdf5-${hdf5_ver}.tar.bz2 ]; then
        echo "hdf5-${hdf5_ver}.tar.bz2 is found"
      else
        download_pkg_from_cp2k_org "${hdf5_sha256}" "hdf5-${hdf5_ver}.tar.bz2"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d hdf5-${hdf5_ver} ] && rm -rf hdf5-${hdf5_ver}
      tar xf hdf5-${hdf5_ver}.tar.bz2
      cd hdf5-${hdf5_ver}
      ./configure \
        --prefix="${pkg_install_dir}" \
        --libdir="${pkg_install_dir}/lib" \
        --enable-fortran \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi
    HDF5_CFLAGS="-I${pkg_install_dir}/include"
    HDF5_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding hdf5 from system paths ===================="
    check_command pkg-config --modversion hdf5
    pkg_install_dir=$(h5cc -show | tr " " "\n" | grep "\-L" | cut -c3-)
    HDF5_CFLAGS="-I${pkg_install_dir}/include"
    HDF5_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking hdf5 to user paths ===================="
    pkg_install_dir="${with_hdf5}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    HDF5_CFLAGS="-I'${pkg_install_dir}/include'"
    HDF5_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "${with_hdf5}" != "__DONTUSE__" ]; then
  # Prefer static libraries if available
  if [ -f "${pkg_install_dir}/lib/libhdf5.a" ]; then
    HDF5_LIBS="-l:libhdf5_fortran.a -l:libhdf5_hl.a -l:libhdf5.a -lz"
  else
    HDF5_LIBS="-lhdf5_fortran -lhdf5_hl -lhdf5 -lz"
  fi
  if [ "${with_hdf5}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_hdf5"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  else
    HDF5_LIBS="${HDF5_LIBS} -lsz"
  fi
  cat << EOF >> "${BUILDDIR}/setup_hdf5"
export HDF5_CFLAGS="${HDF5_CFLAGS}"
export HDF5_LDFLAGS="${HDF5_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__HDF5"
export CP_CFLAGS="\${CP_CFLAGS} ${HDF5_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${HDF5_LDFLAGS}"
export CP_LIBS="${HDF5_LIBS} \${CP_LIBS}"
export HDF5_ROOT="${pkg_install_dir}"
export HDF5_LIBRARIES="${HDF5_LIBS}"
export HDF5_HL_LIBRARIES="${HDF5_LIBS}"
export HDF5_INCLUDE_DIRS="${pkg_install_dir}/include"
EOF
  cat "${BUILDDIR}/setup_hdf5" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_hdf5"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "hdf5"
