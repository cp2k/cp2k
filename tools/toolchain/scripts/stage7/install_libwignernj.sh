#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libwignernj_ver="0.8.0"
libwignernj_sha256="7220cea92652040d6456aba92ff151124d9c69ce8695840490c18dd25a0da80c"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libwignernj" ] && rm "${BUILDDIR}/setup_libwignernj"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libwignernj:=__INSTALL__}" in
  __INSTALL__)
    echo "==================== Installing libwignernj ===================="
    pkg_install_dir="${INSTALLDIR}/libwignernj-${libwignernj_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libwignernj-${libwignernj_ver} is already installed, skipping it."
    else
      retrieve_package "${libwignernj_sha256}" "libwignernj-${libwignernj_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libwignernj-${libwignernj_ver} ] && rm -rf libwignernj-${libwignernj_ver}
      tar -xzf libwignernj-${libwignernj_ver}.tar.gz

      mkdir "libwignernj-${libwignernj_ver}/build"
      cd "libwignernj-${libwignernj_ver}/build"
      # CP2K binds to the C interface through ISO_C_BINDING, so neither the
      # Fortran interface library nor the test suite is needed here
      cmake \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DWIGNERNJ_BUILD_FORTRAN=OFF \
        -DWIGNERNJ_BUILD_TESTS=OFF \
        -DWIGNERNJ_BUILD_CXX_TESTS=OFF \
        -DWIGNERNJ_BUILD_EXAMPLES=OFF \
        .. > cmake.log 2>&1 || tail_excerpt cmake.log
      make -j "$(get_nprocs)" install > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename "${SCRIPT_NAME}")"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding libwignernj from system paths ===================="
    check_lib -lwignernj "libwignernj"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libwignernj.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__)
    report_error "libwignernj is a mandatory dependency of CP2K and cannot be disabled"
    exit 1
    ;;

  *)
    echo "==================== Linking libwignernj to user paths ===================="
    pkg_install_dir="${with_libwignernj}"
    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    LIBWIGNERNJ_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBWIGNERNJ_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${LIBWIGNERNJ_LIBDIR}"
    ;;
esac

if [ "$with_libwignernj" != "__SYSTEM__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libwignernj"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
fi
cat << EOF >> "${BUILDDIR}/setup_libwignernj"
export LIBWIGNERNJ_VER="${libwignernj_ver}"
export LIBWIGNERNJ_ROOT="${pkg_install_dir}"
EOF
filter_setup "${BUILDDIR}/setup_libwignernj" "${SETUPFILE}"

load "${BUILDDIR}/setup_libwignernj"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libwignernj"
