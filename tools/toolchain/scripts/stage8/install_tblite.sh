#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

tblite_ver="0.7.0"
tblite_sha256="3a7cb4602101e828caf41c38ca5e30f82de82d0d26d5db40168acdcad3462b92"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_tblite" ] && rm "${BUILDDIR}/setup_tblite"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_tblite" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing tblite ===================="
    pkg_install_dir="${INSTALLDIR}/tblite-${tblite_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "tblite-${tblite_ver} is already installed, skipping it."
    else
      retrieve_package "${tblite_sha256}" "tblite-${tblite_ver}.tar.xz"
      [ -d tblite-${tblite_ver} ] && rm -rf tblite-${tblite_ver}
      tar -xJf tblite-${tblite_ver}.tar.xz
      cd tblite-${tblite_ver}
      mkdir -p build && cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_TESTING=OFF \
        -DWITH_TESTS=OFF \
        .. \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
      cd ..
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding tblite from system paths ===================="
    check_pkgconfig tblite
    pkg_install_dir="$(pkg-config --variable=prefix tblite)"
    ;;
  *)
    echo "==================== Linking TBLITE to user paths ===================="
    pkg_install_dir="${with_tblite}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac

if [ "$with_tblite" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_tblite"
export TBLITE_ROOT="${pkg_install_dir}"
export TBLITE_VER="${tblite_ver}"
EOF
  if [ "$with_tblite" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_tblite"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_tblite" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_tblite"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "tblite"
