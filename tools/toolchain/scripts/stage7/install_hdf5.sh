#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

hdf5_ver="2.1.1"
hdf5_sha256="efff93b5a904d66e8f626d7da60b5eedc9faf544be27dbabbaa87967b8ad798b"

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
    echo "==================== Installing HDF5 ===================="
    pkg_install_dir="${INSTALLDIR}/hdf5-${hdf5_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "hdf5-${hdf5_ver} is already installed, skipping it."
    else
      retrieve_package "${hdf5_sha256}" "hdf5-${hdf5_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d hdf5-${hdf5_ver} ] && rm -rf hdf5-${hdf5_ver}
      tar xf hdf5-${hdf5_ver}.tar.gz
      cd hdf5-${hdf5_ver}
      mkdir build
      cd build
      # Add "-DHDF5_ENABLE_ZLIB_SUPPORT=ON" because HDF5 2.x doesn't enable ZLIB support by default
      CMAKE_OPTIONS="-DBUILD_TESTING=OFF -DHDF5_BUILD_FORTRAN=ON"
      if [ "$(find_in_paths "libz.*" $LIB_PATHS)" != "__FALSE__" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DHDF5_ENABLE_ZLIB_SUPPORT=ON"
      fi
      if [ "$(find_in_paths "libsz.*" $LIB_PATHS)" != "__FALSE__" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DHDF5_ENABLE_SZIP_SUPPORT=ON"
      fi
      if [ "${MPI_MODE}" != "no" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DHDF5_ENABLE_PARALLEL=ON"
      fi
      cmake .. \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        ${CMAKE_OPTIONS} > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding HDF5 from system paths ===================="
    check_command pkg-config --modversion hdf5
    pkg_install_dir=$(h5cc -showconfig | grep "Installation point" | awk '{print $3}')
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking HDF5 to user paths ===================="
    pkg_install_dir="${with_hdf5}"
    HDF5_LIBDIR="${pkg_install_dir}/lib/MiMiC"
    [ -d "${pkg_install_dir}/lib64" ] && HDF5_LIBDIR="${pkg_install_dir}/lib64/MiMiC"
    check_dir "${HDF5_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "${with_hdf5}" != "__DONTUSE__" ]; then
  if [ "${with_hdf5}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_hdf5"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_hdf5"
export HDF5_VER="${hdf5_ver}"
export HDF5_ROOT="${pkg_install_dir}"
EOF
  filter_setup "${BUILDDIR}/setup_hdf5" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_hdf5"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "hdf5"
