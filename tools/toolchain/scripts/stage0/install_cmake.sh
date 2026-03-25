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

[ -f "${BUILDDIR}/setup_cmake" ] && rm "${BUILDDIR}/setup_cmake"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_cmake}" in
  __INSTALL__)
    echo "==================== Installing CMake ===================="
    cmake_ver="4.3.0"
    if [ "${OPENBLAS_ARCH}" = "arm64" ]; then
      if [ "$(uname -s)" = "Darwin" ]; then
        cmake_arch="macos-universal"
        cmake_sha256="5bd933daf6e9234a53a9a43092746993870d9f162b6c399fd6e4a05cdd475e67"
      elif [ "$(uname -s)" = "Linux" ]; then
        cmake_arch="linux-aarch64"
        cmake_sha256="26fe3011f497eb9398115dcabcc094685e634b1841f7c01dc01c5a89b8b0ea0d"
      else
        report_error ${LINENO} \
          "cmake installation for ARCH=${OPENBLAS_ARCH} under $(uname -s) is not supported. You can try to use the system installation using the flag --with-cmake=system instead."
      fi
    elif [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
      cmake_arch="linux-x86_64"
      cmake_sha256="201bdabe17a54e017f119cffa247648e9c44327e52473c2cc60a88fded94652a"
    else
      report_error ${LINENO} \
        "cmake installation for ARCH=${OPENBLAS_ARCH} under $(uname -s) is not supported. You can try to use the system installation using the flag --with-cmake=system instead."
      exit 1
    fi
    pkg_install_dir="${INSTALLDIR}/cmake-${cmake_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "cmake-${cmake_ver} is already installed, skipping it."
    else
      retrieve_package "${cmake_sha256}" "cmake-${cmake_ver}-${cmake_arch}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      mkdir -p ${pkg_install_dir}
      if [ "${cmake_arch}" = "macos-universal" ]; then
        strip_components=3
      else
        strip_components=1
      fi
      tar --strip-components=${strip_components} -xvf cmake-${cmake_ver}-${cmake_arch}.tar.gz -C ${pkg_install_dir} > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage0/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding CMake from system paths ===================="
    check_command cmake "cmake"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking CMake to user paths ===================="
    pkg_install_dir="${with_cmake}"
    check_dir "${with_cmake}/bin"
    ;;
esac
if [ "${with_cmake}" != "__DONTUSE__" ]; then
  if [ "${with_cmake}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_cmake"
prepend_path PATH "${pkg_install_dir}/bin"
EOF
    filter_setup "${BUILDDIR}/setup_cmake" "${SETUPFILE}"
  fi
fi

load "${BUILDDIR}/setup_cmake"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cmake"
