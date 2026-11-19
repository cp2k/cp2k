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

libint_ver="2.14.1"
libint_pkg="libint-${libint_ver}-cp2k-lmax-${LIBINT_LMAX}.tar.xz"

case "${LIBINT_LMAX}" in
  4)
    libint_sha256="4e71faf2ba1e9b0242c8970c79a74d3b9fa593fdedaaceddb2e2dcd2eed6e314"
    ;;
  5)
    libint_sha256="6e61df7a3bbc298ded738698cb735f01bf43a2d6ec996360c6cd4b31346c571a"
    ;;
  6)
    libint_sha256="72bda0d0eb92f63600c3fdf65806ed9b4d49a0a3c706a3576e9a2c0c2619f31c"
    ;;
  7)
    libint_sha256="aa78b4868050e3e8605e1cb3b96f5a81238ae4b24d83eb01e9845f95c8999eec"
    ;;
  *)
    report_error "Unsupported value --libint-lmax=${LIBINT_LMAX}."
    exit 1
    ;;
esac

[ -f "${BUILDDIR}/setup_libint" ] && rm "${BUILDDIR}/setup_libint"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libint" in
  __INSTALL__)
    echo "==================== Installing LIBINT ===================="
    pkg_install_dir="${INSTALLDIR}/libint-${libint_ver}-cp2k-lmax-${LIBINT_LMAX}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libint-${libint_ver} is already installed, skipping it."
    else
      retrieve_package "${libint_sha256}" "${libint_pkg}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libint-${libint_ver}-cp2k-lmax-${LIBINT_LMAX} ] && rm -rf libint-${libint_ver}-cp2k-lmax-${LIBINT_LMAX}
      tar -xJf ${libint_pkg}
      cd libint-${libint_ver}-cp2k-lmax-${LIBINT_LMAX}
      # reduce debug information to level 1 since
      # level 2 (default for -g flag) leads to very large binary size
      cmake -Bbuild \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DCMAKE_CXX_FLAGS="-g1" \
        -DBUILD_TESTING=OFF \
        -DLIBINT2_ENABLE_FORTRAN=OFF \
        -DLIBINT2_REQUIRE_CXX_API=ON \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      cmake --build build -t install -j $(get_nprocs) > build.log 2>&1 || tail_excerpt build.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding LIBINT from system paths ===================="
    check_lib -lint2 "libint"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libint2.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LIBINT to user paths ===================="
    pkg_install_dir="${with_libint}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_libint" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libint"
export LIBINT2_ROOT="${pkg_install_dir}"
export LIBINT_VER="${libint_ver}"
EOF
  if [ "$with_libint" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libint"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libint" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libint"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libint"
