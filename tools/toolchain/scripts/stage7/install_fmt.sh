#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

fmt_ver="12.1.0"
fmt_sha256="695fd197fa5aff8fc67b5f2bbc110490a875cdf7a41686ac8512fb480fa8ada7"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_fmt" ] && rm "${BUILDDIR}/setup_fmt"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_fmt" in
  __INSTALL__)
    echo "==================== Installing fmt ===================="
    pkg_install_dir="${INSTALLDIR}/fmt-${fmt_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "fmt-${fmt_ver} is already installed, skipping it."
    else
      retrieve_package "${fmt_sha256}" "fmt-${fmt_ver}.zip"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d fmt-${fmt_ver} ] && rm -rf fmt-${fmt_ver}
      unzip -q fmt-${fmt_ver}.zip
      cd fmt-${fmt_ver}
      mkdir build
      cd build
      cmake .. \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DFMT_TEST=OFF \
        > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding fmt from system paths ===================="
    check_lib -lfmt "fmt"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking fmt to user paths ===================="
    pkg_install_dir="${with_fmt}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "${with_fmt}" != "__DONTUSE__" ]; then
  if [ "${with_fmt}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_fmt"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_fmt"
export fmt_VER="${fmt_ver}"
export fmt_ROOT="${pkg_install_dir}"
EOF
  filter_setup "${BUILDDIR}/setup_fmt" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_fmt"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "fmt"
