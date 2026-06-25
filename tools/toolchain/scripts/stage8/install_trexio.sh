#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

trexio_ver="2.6.1"
trexio_sha256="c3694ec1528632a386a2af89199c75d70ecd45bfcc2ca1d4ccccbfa1308ad5fa"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_trexio" ] && rm "${BUILDDIR}/setup_trexio"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_trexio" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing TREXIO ===================="
    require_env HDF5_ROOT

    pkg_install_dir="${INSTALLDIR}/trexio-${trexio_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "trexio-${trexio_ver} is already installed, skipping it."
    else
      retrieve_package "${trexio_sha256}" "trexio-${trexio_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d trexio-${trexio_ver} ] && rm -rf trexio-${trexio_ver}
      tar -xzf trexio-${trexio_ver}.tar.gz
      cd trexio-${trexio_ver}

      CMAKE_OPTIONS="-DCMAKE_VERBOSE_MAKEFILE=ON"
      if [ "${MPI_MODE}" != "no" ]; then
        CMAKE_OPTIONS+="-DCMAKE_C_COMPILER=${MPICC}"
      fi

      mkdir build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        ${CMAKE_OPTIONS} \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) >> make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding trexio from system paths ===================="
    require_env HDF5_ROOT
    check_lib -ltrexio "trexio"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libtrexio.*" $LIB_PATHS)))"
    ;;
  *)
    echo "==================== Linking TREXIO to user paths ===================="
    pkg_install_dir="${with_trexio}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_trexio" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_trexio"
export TREXIO_VER="${trexio_ver}"
export TREXIO_ROOT="${pkg_install_dir}"
EOF
  if [ "$with_trexio" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_trexio"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_trexio" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_trexio"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "trexio"
