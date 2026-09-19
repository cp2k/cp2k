#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

libwignernj_ver="0.8.0"
libwignernj_sha256="7220cea92652040d6456aba92ff151124d9c69ce8695840490c18dd25a0da80c"

source "${SCRIPT_DIR}/common_vars.sh"
source "${SCRIPT_DIR}/tool_kit.sh"
source "${SCRIPT_DIR}/signal_trap.sh"
source "${INSTALLDIR}/toolchain.conf"
source "${INSTALLDIR}/toolchain.env"

rm -f "${BUILDDIR}/setup_libwignernj"
mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libwignernj:=__INSTALL__}" in
  __INSTALL__)
    echo "==================== Installing libwignernj ===================="
    pkg_install_dir="${INSTALLDIR}/libwignernj-${libwignernj_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libwignernj-${libwignernj_ver} is already installed, skipping it."
    else
      archive="libwignernj-${libwignernj_ver}.tar.gz"
      if ! [ -f "${archive}" ] || ! checksum "${libwignernj_sha256}" "${archive}"; then
        download_pkg_from_urlpath "${libwignernj_sha256}" "v${libwignernj_ver}.tar.gz" \
          "https://github.com/susilehtola/libwignernj/archive/refs/tags" "${archive}"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf "libwignernj-${libwignernj_ver}"
      tar -xzf "${archive}"
      mkdir "libwignernj-${libwignernj_ver}/build"
      cd "libwignernj-${libwignernj_ver}/build"
      cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DWIGNERNJ_BUILD_FORTRAN=OFF \
        -DWIGNERNJ_BUILD_TESTS=OFF \
        -DWIGNERNJ_BUILD_CXX_TESTS=OFF \
        -DWIGNERNJ_BUILD_EXAMPLES=OFF \
        -DWIGNERNJ_BUILD_LTO=OFF \
        .. > cmake.log 2>&1 || tail_excerpt cmake.log
      make -j "$(get_nprocs)" install > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/install_libwignernj.sh"
    fi
    LIBWIGNERNJ_LIBDIR="${pkg_install_dir}/lib"
    ;;
  __SYSTEM__)
    echo "==================== Finding libwignernj from system paths ===================="
    check_lib -lwignernj "libwignernj"
    library_path=$(find_in_paths "libwignernj.*" ${LIB_PATHS})
    if [ "${library_path}" = "__FALSE__" ]; then
      report_error "Cannot locate libwignernj in the system library paths"
    fi
    LIBWIGNERNJ_LIBDIR=$(dirname "${library_path}")
    pkg_install_dir=$(dirname "${LIBWIGNERNJ_LIBDIR}")
    ;;
  __DONTUSE__)
    report_error "libwignernj is a required dependency of CP2K and cannot be disabled"
    ;;
  *)
    echo "==================== Linking libwignernj to user paths ===================="
    pkg_install_dir="${with_libwignernj}"
    LIBWIGNERNJ_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBWIGNERNJ_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${LIBWIGNERNJ_LIBDIR}"
    ;;
esac

if [ "${with_libwignernj}" != "__SYSTEM__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libwignernj"
prepend_path LD_LIBRARY_PATH "${LIBWIGNERNJ_LIBDIR}"
prepend_path LD_RUN_PATH "${LIBWIGNERNJ_LIBDIR}"
prepend_path LIBRARY_PATH "${LIBWIGNERNJ_LIBDIR}"
prepend_path PKG_CONFIG_PATH "${LIBWIGNERNJ_LIBDIR}/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
fi

cat << EOF >> "${BUILDDIR}/setup_libwignernj"
export wignernj_ROOT="${pkg_install_dir}"
EOF
filter_setup "${BUILDDIR}/setup_libwignernj" "${SETUPFILE}"
load "${BUILDDIR}/setup_libwignernj"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libwignernj"
