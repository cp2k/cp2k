#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

skala_model_ver="1.1"
skala_model_rev="rev1"
skala_model_pkg="skala-${skala_model_ver}-${skala_model_rev}.fun"
skala_model_sha256="7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd"
skala_model_urlpath="https://huggingface.co/microsoft/skala-${skala_model_ver}/resolve/main"
skala_cuda_model_pkg="skala-${skala_model_ver}-${skala_model_rev}-cuda.fun"
skala_cuda_model_sha256="f848eae769dca91741a518ae7275d10caac398ab21db649f91bc1f136872f223"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_skala" ] && rm "${BUILDDIR}/setup_skala"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_skala}" in
  __INSTALL__)
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
      # Download and verify non-CUDA model
      if ! [ -f "${skala_model_pkg}" ]; then
        echo "Downloading skala model..."
        download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
      fi
      if ! checksum "${skala_model_sha256}" "${skala_model_pkg}"; then
        echo "${skala_model_pkg} checksum failed, re-downloading..."
        rm -f "${skala_model_pkg}"
        download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
      fi

      # Download and verify CUDA model if ENABLE_CUDA
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
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage5/$(basename "${SCRIPT_NAME}")" \
          "${skala_model_pkg}" "${skala_cuda_model_pkg}"
      else
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage5/$(basename "${SCRIPT_NAME}")" \
          "${skala_model_pkg}"
      fi
    fi

    # Create setup file (whether newly installed or already existing)
    cat > "${BUILDDIR}/setup_skala" << EOF
SKALA_MODEL_VER="${skala_model_ver}"
SKALA_MODEL_ROOT="${skala_install_dir}"
SKALA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
EOF
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      cat >> "${BUILDDIR}/setup_skala" << EOF
SKALA_CUDA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
EOF
    fi

    filter_setup "${BUILDDIR}/setup_skala" "${SETUPFILE}"
    load "${BUILDDIR}/setup_skala"
    write_toolchain_env "${INSTALLDIR}"
    ;;
  __SYSTEM__)
    # Skala is expected to be provided by system, inform user about required environment variables
    echo "Skala installation skipped - expecting system installation (with_skala=${with_skala})"
    echo "The following environment variables should be set for system Skala:"
    echo "  SKALA_MODEL_ROOT - Path to Skala installation directory"
    echo "  SKALA_MODEL - Path to non-CUDA Skala model file"
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      echo "  SKALA_CUDA_MODEL - Path to CUDA Skala model file (required when CUDA is enabled)"
    fi
    ;;
  __DONTUSE__)
    # Skala is explicitly disabled, do nothing
    ;;
  *)
    # Custom Skala installation path
    echo "==================== Linking Skala to custom path ===================="
    skala_install_dir="${with_skala}"
    check_dir "${skala_install_dir}/share/skala/onedft_models"
    check_dir "${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      check_dir "${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
    fi
    # Create setup file for custom path
    cat > "${BUILDDIR}/setup_skala" << EOF
SKALA_MODEL_VER="${skala_model_ver}"
SKALA_MODEL_ROOT="${skala_install_dir}"
SKALA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
EOF
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      cat >> "${BUILDDIR}/setup_skala" << EOF
SKALA_CUDA_MODEL="${skala_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
EOF
    fi
    filter_setup "${BUILDDIR}/setup_skala" "${SETUPFILE}"
    load "${BUILDDIR}/setup_skala"
    ;;
esac

# For __SYSTEM__ and __DONTUSE__ cases, create appropriate setup_skala file
case "${with_skala}" in
  __SYSTEM__)
    # For system installation, create setup file with expected environment variables
    cat > "${BUILDDIR}/setup_skala" << EOF
SKALA_MODEL_VER="${skala_model_ver}"
SKALA_MODEL_ROOT="${SKALA_MODEL_ROOT:-}"
SKALA_MODEL="${SKALA_MODEL:-}"
EOF
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      cat >> "${BUILDDIR}/setup_skala" << EOF
SKALA_CUDA_MODEL="${SKALA_CUDA_MODEL:-}"
EOF
    fi
    filter_setup "${BUILDDIR}/setup_skala" "${SETUPFILE}"
    load "${BUILDDIR}/setup_skala"
    ;;
  __DONTUSE__)
    # For disabled skala, create empty setup file
    cat > "${BUILDDIR}/setup_skala" << EOF
SKALA_MODEL_VER=""
SKALA_MODEL_ROOT=""
SKALA_MODEL=""
EOF
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      cat >> "${BUILDDIR}/setup_skala" << EOF
SKALA_CUDA_MODEL=""
EOF
    fi
    filter_setup "${BUILDDIR}/setup_skala" "${SETUPFILE}"
    load "${BUILDDIR}/setup_skala"
    ;;
  *)
    # __INSTALL__ and custom path cases already handled above
    ;;
esac

cd "${ROOTDIR}"
report_timing "skala"
