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

eigen_ver="5.0.1"
eigen_sha256="e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec"

[ -f "${BUILDDIR}/setup_eigen" ] && rm "${BUILDDIR}/setup_eigen"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_eigen" in
  __INSTALL__)
    echo "==================== Installing Eigen ===================="
    pkg_install_dir="${INSTALLDIR}/eigen-${eigen_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "eigen-${eigen_ver} is already installed, skipping it."
    else
      retrieve_package "${eigen_sha256}" "eigen-${eigen_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d eigen-${eigen_ver} ] && rm -rf eigen-${eigen_ver}
      tar -xzf "eigen-${eigen_ver}.tar.gz"
      cd eigen-${eigen_ver}

      mkdir build
      cd build
      cmake .. \
        -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        -DBUILD_TESTING=OFF \
        -DEIGEN_BUILD_BLAS=OFF \
        -DEIGEN_BUILD_LAPACK=OFF \
        > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log

      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Eigen from system paths ===================="
    check_pkgconfig eigen3
    pkg_install_dir="$(pkg-config --variable=prefix eigen3)"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Eigen to user paths ===================="
    pkg_install_dir="$with_eigen"
    check_dir "${pkg_install_dir}/include/eigen3"
    ;;
esac
if [ "$with_eigen" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_eigen"
export EIGEN_VER="${eigen_ver}"
export EIGEN_ROOT="${pkg_install_dir}"
EOF
  if [ "$with_eigen" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_eigen"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
prepend_path CPATH "${pkg_install_dir}/include/eigen3"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_eigen" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_eigen"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "eigen"
