#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# From https://pytorch.org/get-started/locally/
libtorch_ver="2.4.0"
libtorch_sha256="63d572598c8d532128a335018913e795c1bbb32602ce378896dc8cfbb5590976"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libtorch" ] && rm "${BUILDDIR}/setup_libtorch"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libtorch}" in
  __INSTALL__)
    echo "==================== Installing libtorch ===================="
    pkg_install_dir="${INSTALLDIR}/libtorch-${libtorch_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    
    # Determine CUDA suffix based on ENABLE_CUDA
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      # Use CUDA 12.1 for A100/A40/H100, CUDA 11.8 for older GPUs
      if [ "${GPUVER}" = "A100" ] || [ "${GPUVER}" = "A40" ] || [ "${GPUVER}" = "H100" ]; then
        cuda_suffix="+cu121"
        # SHA256 for libtorch-cxx11-abi-shared-with-deps-2.4.0+cu121.zip
        # From PyTorch official downloads (check for updates)
        libtorch_sha256="8f1f7063e421bcc627f9e2b0779331eddb8c8417ce657ca8c25c9f8b6b3cf2a0"
      else
        cuda_suffix="+cu118"
        # SHA256 for libtorch-cxx11-abi-shared-with-deps-2.4.0+cu118.zip
        # From PyTorch official downloads (check for updates)
        libtorch_sha256="f739db778882e8826b92ab9e140c9c66a05041c621121386aae718c0110679fc"
      fi
      archive_file="libtorch-cxx11-abi-shared-with-deps-${libtorch_ver}${cuda_suffix}.zip"
      echo "Installing CUDA-enabled libtorch (${GPUVER})"
    else
      archive_file="libtorch-cxx11-abi-shared-with-deps-${libtorch_ver}+cpu.zip"
      echo "Installing CPU-only libtorch"
    fi

    if verify_checksums "${install_lock_file}"; then
      echo "libtorch-${libtorch_ver} is already installed, skipping it."
    else
      # Use PyTorch official download URL for CUDA builds
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        if [ "${GPUVER}" = "A100" ] || [ "${GPUVER}" = "A40" ] || [ "${GPUVER}" = "H100" ]; then
          download_url="https://download.pytorch.org/libtorch/cu121/${archive_file}"
        else
          download_url="https://download.pytorch.org/libtorch/cu118/${archive_file}"
        fi
        echo "Downloading from: ${download_url}"
        wget --quiet "${download_url}" -O "${archive_file}" || \
          report_error ${LINENO} "Failed to download ${archive_file} from ${download_url}"
      else
        retrieve_package "${libtorch_sha256}" "${archive_file}"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libtorch ] && rm -rf libtorch
      [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
      unzip -q ${archive_file}
      mv libtorch ${pkg_install_dir}

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename "${SCRIPT_NAME}")"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding libtorch from system paths ===================="
    check_lib -ltorch "libtorch"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libtorch.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking libtorch to user paths ===================="
    pkg_install_dir="${with_libtorch}"
    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    LIBTORCH_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBTORCH_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${LIBTORCH_LIBDIR}"
    ;;
esac

if [ "$with_libtorch" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libtorch"
export LIBTORCH_VER="${libtorch_ver}"
EOF
  if [ "$with_libtorch" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libtorch"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libtorch" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libtorch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libtorch"
