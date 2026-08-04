#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")" && pwd -P)"

skala_model_ver="1.1"
skala_model_pkg="skala-${skala_model_ver}.fun"
skala_model_sha256="0c8432ac3f03c8f1276372df9aca5b7ee7f8939d47a8789eb158976e89aa0606"
skala_model_urlpath="https://huggingface.co/microsoft/skala-${skala_model_ver}/resolve/main"

source "${SCRIPT_DIR}/.."/common_vars.sh
source "${SCRIPT_DIR}/.."/tool_kit.sh
source "${SCRIPT_DIR}/.."/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_skala" ] && rm "${BUILDDIR}/setup_skala"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

# Shared installation directory for skala model
skala_install_dir="${INSTALLDIR}/skala-${skala_model_ver}"
install_lock_file="${skala_install_dir}/install_successful"

echo "==================== Installing Skala model ===================="

# Check if already installed
if verify_checksums "${install_lock_file}"; then
  echo "skala-${skala_model_ver} is already installed, skipping it."
else
  # Download skala model if not already present
  if ! [ -f "${skala_model_pkg}" ]; then
    echo "Downloading skala model..."
    download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
  fi

  # Verify checksum
  if ! checksum "${skala_model_sha256}" "${skala_model_pkg}"; then
    echo "${skala_model_pkg} checksum failed, re-downloading..."
    rm -f "${skala_model_pkg}"
    download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
  fi

  # Create installation directory
  mkdir -p "${skala_install_dir}/share/skala/onedft_models"

  # Install skala model
  if ! [ -f "${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}" ]; then
    echo "Installing skala model to ${skala_install_dir}..."
    install -m 0644 "${skala_model_pkg}" "${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
  fi

  # Create setup_skala file
  cat > "${BUILDDIR}/setup_skala" << EOF
SKALA_MODEL_VER="${skala_model_ver}"
SKALA_MODEL_ROOT="${skala_install_dir}"
SKALA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
EOF

  echo "Skala model installed successfully to ${skala_install_dir}"

  # Write checksums for verification
  write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename "${SCRIPT_NAME}")" \
    "${skala_model_pkg}"
fi

# Create setup file (whether newly installed or already existing)
cat > "${BUILDDIR}/setup_skala" << EOF
SKALA_MODEL_VER="${skala_model_ver}"
SKALA_MODEL_ROOT="${skala_install_dir}"
SKALA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
EOF

filter_setup "${BUILDDIR}/setup_skala" "${SETUPFILE}"
load "${BUILDDIR}/setup_skala"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "skala"
