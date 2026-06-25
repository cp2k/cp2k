#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

deepmd_ver="3.1.0"
deepmd_sha256="45f13df9ed011438d139a7f61416b8d7940f63c47fcde53180bfccd60c9d22ee"
deepmd_pkg="deepmd-kit-${deepmd_ver}.tar.gz"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_deepmd" ] && rm "${BUILDDIR}/setup_deepmd"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_deepmd" in
  __INSTALL__)
    echo "==================== Installing DeePMD ===================="
    pkg_install_dir="${INSTALLDIR}/deepmd-kit-${deepmd_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libdeepmd_c-${deepmd_ver} is already installed, skipping it."
    else
      retrieve_package "${deepmd_sha256}" "${deepmd_pkg}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d deepmd-kit-${deepmd_ver} ] && rm -rf deepmd-kit-${deepmd_ver}
      tar -xzf ${deepmd_pkg}
      cd deepmd-kit-${deepmd_ver}/source

      mkdir build
      cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_CXX_STANDARD=11 \
        -DCMAKE_CXX_STANDARD_REQUIRED=TRUE \
        -DENABLE_PYTORCH=TRUE \
        .. > cmake.log 2>&1 || tail_excerpt cmake.log
      make -j deepmd_c > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding DeePMD from system paths ===================="
    check_lib -ldeepmd "DEEPMD"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libdeepmd.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking DEEPMD to user paths ===================="
    pkg_install_dir="${with_deepmd}"
    check_dir "${pkg_install_dir}/include/deepmd"
    check_dir "${pkg_install_dir}/lib"
    ;;
esac

if [ "$with_deepmd" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_deepmd"
export DEEPMD_VER="${deepmd_ver}"
export DEEPMD_ROOT="${pkg_install_dir}"
EOF
  if [ "$with_deepmd" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_deepmd"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
    filter_setup "${BUILDDIR}/setup_deepmd" "${SETUPFILE}"
  fi
  cat << EOF >> "${INSTALLDIR}/lsan.supp"
# leaks related to DeePMD-kit and TensorFlow
leak:DP_NewDeepPot
leak:deepmd::AtomMap::AtomMap
EOF
fi

load "${BUILDDIR}/setup_deepmd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "deepmd"
