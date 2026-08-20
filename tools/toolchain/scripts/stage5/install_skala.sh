#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")" && pwd -P)"

skala_model_ver="1.1"
skala_model_rev="rev1"
skala_model_pkg="skala-${skala_model_ver}-${skala_model_rev}.fun"
skala_model_sha256="7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd"
skala_model_urlpath="https://huggingface.co/microsoft/skala-${skala_model_ver}/resolve/main"
skala_cuda_model_pkg="skala-${skala_model_ver}-${skala_model_rev}-cuda.fun"
skala_cuda_model_sha256="f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223"

source "${SCRIPT_DIR}/.."/common_vars.sh
source "${SCRIPT_DIR}/.."/tool_kit.sh
source "${SCRIPT_DIR}/.."/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${INSTALLDIR}/setup_skala" ] && rm "${INSTALLDIR}/setup_skala"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

# Shared installation directory for skala model
skala_install_dir="${INSTALLDIR}/skala-${skala_model_ver}"
install_lock_file="${skala_install_dir}/install_successful"

echo "==================== Installing Skala model ===================="

# Check if already installed: non-CUDA model always needed; CUDA model only if ENABLE_CUDA
non_cuda_ok() {
  [ -f "${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}" ]
}
cuda_ok() {
  [ "${ENABLE_CUDA}" != "__TRUE__" ] ||
    [ -f "${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}" ]
}

if non_cuda_ok && cuda_ok && verify_checksums "${install_lock_file}"; then
  echo "skala-${skala_model_ver} is already installed, skipping it."
else
  # Download non-CUDA model
  if ! [ -f "${skala_model_pkg}" ]; then
    echo "Downloading skala model..."
    download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
  fi
  if ! checksum "${skala_model_sha256}" "${skala_model_pkg}"; then
    echo "${skala_model_pkg} checksum failed, re-downloading..."
    rm -f "${skala_model_pkg}"
    download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
  fi

  # Download CUDA model if ENABLE_CUDA
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    if ! [ -f "${skala_cuda_model_pkg}" ]; then
      echo "Downloading skala CUDA model..."
      download_pkg_from_urlpath "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}" "${skala_model_urlpath}" "${skala_cuda_model_pkg}"
    fi
    if ! checksum "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}"; then
      echo "${skala_cuda_model_pkg} checksum failed, re-downloading..."
      rm -f "${skala_cuda_model_pkg}"
      download_pkg_from_urlpath "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}" "${skala_model_urlpath}" "${skala_cuda_model_pkg}"
    fi
  fi

  # Create installation directory
  mkdir -p "${skala_install_dir}/share/skala/onedft_models"

  # Install non-CUDA model
  if ! [ -f "${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}" ]; then
    echo "Installing skala model to ${skala_install_dir}..."
    install -m 0644 "${skala_model_pkg}" "${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
  fi

  # Install CUDA model if ENABLE_CUDA
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    if ! [ -f "${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}" ]; then
      echo "Installing skala CUDA model to ${skala_install_dir}..."
      install -m 0644 "${skala_cuda_model_pkg}" "${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
    fi
  fi

  echo "Skala model installed successfully to ${skala_install_dir}"

  # Write checksums for verification
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename "${SCRIPT_NAME}")" \
      "${skala_model_pkg}" "${skala_cuda_model_pkg}"
  else
    write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename "${SCRIPT_NAME}")" \
      "${skala_model_pkg}"
  fi
fi

# Create setup file (whether newly installed or already existing)
cat > "${INSTALLDIR}/setup_skala" << EOF
SKALA_MODEL_VER="${skala_model_ver}"
SKALA_MODEL_ROOT="${skala_install_dir}"
SKALA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
EOF
if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
  cat >> "${INSTALLDIR}/setup_skala" << EOF
SKALA_CUDA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
EOF
fi

filter_setup "${INSTALLDIR}/setup_skala" "${SETUPFILE}"
load "${INSTALLDIR}/setup_skala"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "skala"
