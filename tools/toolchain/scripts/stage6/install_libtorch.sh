#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# From https://pytorch.org/get-started/locally/
libtorch_ver="2.7.1"
libtorch_sha256="63d572598c8d532128a335018913e795c1bbb32602ce378896dc8cfbb5590976"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libtorch" ] && rm "${BUILDDIR}/setup_libtorch"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libtorch}" in
  __INSTALL__)
    echo "==================== Installing libtorch ===================="
    pkg_install_dir="${INSTALLDIR}/libtorch-${libtorch_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    archive_file="libtorch-cxx11-abi-shared-with-deps-${libtorch_ver}+cpu.zip"

    if verify_checksums "${install_lock_file}"; then
      echo "libtorch-${libtorch_ver} is already installed, skipping it."
    else
      retrieve_package "${libtorch_sha256}" "${archive_file}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libtorch ] && rm -rf libtorch
      [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
      unzip -q ${archive_file}
      mv libtorch ${pkg_install_dir}

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename "${SCRIPT_NAME}")"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding libtorch from system paths ===================="
    check_lib -ltorch "libtorch"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libtorch.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking libtorch to user paths ===================="
    pkg_install_dir="${with_libtorch}"
    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    LIBTORCH_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBTORCH_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${LIBTORCH_LIBDIR}"
    ;;
esac

if [ "$with_libtorch" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libtorch"
export LIBTORCH_VER="${libtorch_ver}"
EOF
  if [ "$with_libtorch" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libtorch"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libtorch" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libtorch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libtorch"
