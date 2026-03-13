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
    cmake_ver="4.2.3"
    cmake_ext="sh"
    if [ "${OPENBLAS_ARCH}" = "arm64" ]; then
      if [ "$(uname -s)" = "Darwin" ]; then
        cmake_arch="macos-universal"
        cmake_sha256="c2302d3e9c48daabee5ea7c4db4b2b93b989bcc89dae8b760880e00120641b5b"
        cmake_ext="tar.gz"
      elif [ "$(uname -s)" = "Linux" ]; then
        cmake_arch="linux-aarch64"
        cmake_sha256="8e65cb924c41e3f6c5c98d1900aeeab3cfe3494f18ed8e5e63b91022df710401"
      else
        report_error ${LINENO} \
          "cmake installation for ARCH=${OPENBLAS_ARCH} under $(uname -s) is not supported. You can try to use the system installation using the flag --with-cmake=system instead."
      fi
    elif [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
      cmake_arch="linux-x86_64"
      cmake_sha256="b760514fde7fc510fcd16e51a81a4d2687b1f051b263d40b6806789d3d9fd62c"
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
      # retrieve_package "${cmake_sha256}" "cmake-${cmake_ver}-${cmake_arch}.${cmake_ext}"
      if [ -f cmake-${cmake_ver}-${cmake_arch}.${cmake_ext} ]; then
        echo "cmake-${cmake_ver}-${cmake_arch}.${cmake_ext} is found"
      else
        download_pkg_from_urlpath "${cmake_sha256}" "cmake-${cmake_ver}-${cmake_arch}.${cmake_ext}" https://github.com/Kitware/CMake/releases/download/v4.2.3
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      mkdir -p ${pkg_install_dir}
      if [ "${cmake_arch}" = "macos-universal" ]; then
        tar --strip-components=3 -xvf cmake-${cmake_ver}-${cmake_arch}.${cmake_ext} -C ${pkg_install_dir} > install.log 2>&1 || tail_excerpt install.log
      else
        /bin/sh cmake-${cmake_ver}-${cmake_arch}.${cmake_ext} --prefix=${pkg_install_dir} --skip-license > install.log 2>&1 || tail_excerpt install.log
      fi
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
