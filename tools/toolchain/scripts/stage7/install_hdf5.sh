#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

hdf5_ver="1.12.0"
hdf5_sha256="97906268640a6e9ce0cde703d5a71c9ac3092eded729591279bf2e3ca9765f61"

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
        download_pkg ${DOWNLOADER_FLAGS} ${hdf5_sha256} \
          https://www.cp2k.org/static/downloads/hdf5-${hdf5_ver}.tar.bz2
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d hdf5-${hdf5_ver} ] && rm -rf hdf5-${hdf5_ver}
      tar xf hdf5-${hdf5_ver}.tar.bz2
      cd hdf5-${hdf5_ver}
      ./configure \
        --prefix="${pkg_install_dir}" \
        --libdir="${pkg_install_dir}/lib" \
        --enable-shared \
        > configure.log 2>&1
      make -j $(get_nprocs) > make.log 2>&1
      make -j $(get_nprocs) install > install.log 2>&1
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi

    HDF5_CFLAGS="-I${pkg_install_dir}/include"
    HDF5_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding hdf5 from system paths ===================="
    check_command pkg-config --modversion hdf5
    add_include_from_paths HDF5_CFLAGS "hdf5.h" $INCLUDE_PATHS
    add_lib_from_paths HDF5_LDFLAGS "libhdf5.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *) ;;

esac
if [ "$with_hdf5" != "__DONTUSE__" ]; then
  HDF5_LIBS="-lhdf5 -lhdf5_hl"
  if [ "$with_hdf5" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_hdf5"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_hdf5"
export HDF5_CFLAGS="${HDF5_CFLAGS}"
export HDF5_LDFLAGS="${HDF5_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__HDF5|)"
export CP_CFLAGS="\${CP_CFLAGS} ${HDF5_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${HDF5_LDFLAGS}"
####################################################
#
# include hdf5 only if sirius is activated and build
# depends them on mpi
#
####################################################

export CP_LIBS="IF_MPI(${HDF5_LIBS}|) \${CP_LIBS}"
export HDF5_ROOT="$pkg_install_dir"
export HDF5_LIBRARIES="$HDF5_LIBS"
export HDF5_HL_LIBRARIES="$HDF5_LIBS"
export HDF5_INCLUDE_DIRS="$pkg_install_dir/include"

EOF
  cat "${BUILDDIR}/setup_hdf5" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_hdf5"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "hdf5"
