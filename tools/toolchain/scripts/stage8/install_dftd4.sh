#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

dftd4_ver="4.2.0"
dftd4_sha256="467e024071510ad82b862c66c383c2ebc164fc1140e15dfc79f48d2f999fd184"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_dftd4" ] && rm "${BUILDDIR}/setup_dftd4"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_dftd4" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing DFTD4 ===================="
    pkg_install_dir="${INSTALLDIR}/dftd4-${dftd4_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "dftd4-${dftd4_ver} is already installed, skipping it."
    else
      retrieve_package "${dftd4_sha256}" dftd4-${dftd4_ver}.tar.xz
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d dftd4-${dftd4_ver} ] && rm -rf dftd4-${dftd4_ver}
      tar -xJf dftd4-${dftd4_ver}.tar.xz
      cd dftd4-${dftd4_ver}

      mkdir build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        .. \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding DFTD4 from system paths ===================="
    check_command pkg-config --modversion dftd4
    pkg_install_dir="$(pkg-config --variable=prefix dftd4)"
    ;;
  *)
    echo "==================== Linking DFTD4 to user paths ===================="
    pkg_install_dir="$with_dftd4"
    check_dir "${pkg_install_dir}/include"
    ;;
esac

if [ "$with_dftd4" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_dftd4"
export DFTD4_ROOT="${pkg_install_dir}"
export DFTD4_VER="${dftd4_ver}"
EOF
  if [ "$with_dftd4" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_dftd4"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_dftd4" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_dftd4"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "dftd4"
