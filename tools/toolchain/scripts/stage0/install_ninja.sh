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

[ -f "${BUILDDIR}/setup_ninja" ] && rm "${BUILDDIR}/setup_ninja"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_ninja}" in
  __INSTALL__)
    echo "==================== Installing Ninja  ===================="
    ninja_ver="1.12.1"
    ninja_sha256="821bdff48a3f683bc4bb3b6f0b5fe7b2d647cf65d52aeb63328c91a6c6df285a"

    pkg_install_dir="${INSTALLDIR}/ninja-v${ninja_ver}"
    install_lock_file="$pkg_install_dir/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "ninja-v${ninja_ver} is already installed, skipping it."
    else
      if [ -f ninja-v${ninja_ver}.tar.gz ]; then
        echo "ninja-v${ninja_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${ninja_sha256}" "ninja-v${ninja_ver}.tar.gz"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      mkdir -p ${pkg_install_dir}
      tar -xzf ninja-v${ninja_ver}.tar.gz
      cd ninja-${ninja_ver}
      cmake -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        -Bbuild-ninja \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      cmake --build build-ninja -j $(get_nprocs) > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      cmake --install build-ninja > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage0/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Ninja from system paths ===================="
    check_command ninja "ninja"
    ;;
  __DONTUSE__)
    # Nothing to do
    echo "Ninja required for DFTD4"
    ;;
  *)
    echo "==================== Linking Ninja to user paths ===================="
    pkg_install_dir="${with_ninja}"
    check_dir "${with_ninja}/bin"
    ;;
esac
if [ "${with_ninja}" != "__DONTUSE__" ]; then
  if [ "${with_ninja}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_ninja"
prepend_path PATH "${pkg_install_dir}/bin"
EOF
    cat "${BUILDDIR}/setup_ninja" >> $SETUPFILE
  fi
fi

load "${BUILDDIR}/setup_ninja"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "ninja"
