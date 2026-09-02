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
    echo "==================== Installing Skala model ===================="
    pkg_install_dir="${INSTALLDIR}/skala-${skala_model_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    install_ok() {
      [ -f "${pkg_install_dir}/share/skala/onedft_models/${skala_model_pkg}" ] &&
        { [ "${ENABLE_CUDA}" != "__TRUE__" ] ||
          [ -f "${pkg_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}" ]; }
    }

    if install_ok && verify_checksums "${install_lock_file}"; then
      echo "skala-${skala_model_ver} is already installed, skipping it."
    else
      mkdir -p "${pkg_install_dir}/share/skala/onedft_models"

      echo "Downloading skala model..."
      if [ -f "${skala_model_pkg}" ] && checksum "${skala_model_sha256}" "${skala_model_pkg}"; then
        :
      else
        download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
        if ! checksum "${skala_model_sha256}" "${skala_model_pkg}"; then
          echo "${skala_model_pkg} checksum failed, re-downloading..."
          rm -f "${skala_model_pkg}"
          download_pkg_from_urlpath "${skala_model_sha256}" "${skala_model_pkg}" "${skala_model_urlpath}" "${skala_model_pkg}"
        fi
      fi
      echo "Installing skala model to ${pkg_install_dir}..."
      install -m 0644 "${skala_model_pkg}" "${pkg_install_dir}/share/skala/onedft_models/${skala_model_pkg}"

      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        echo "Downloading skala CUDA model..."
        if [ -f "${skala_cuda_model_pkg}" ] && checksum "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}"; then
          :
        else
          download_pkg_from_urlpath "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}" "${skala_model_urlpath}" "${skala_cuda_model_pkg}"
          if ! checksum "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}"; then
            echo "${skala_cuda_model_pkg} checksum failed, re-downloading..."
            rm -f "${skala_cuda_model_pkg}"
            download_pkg_from_urlpath "${skala_cuda_model_sha256}" "${skala_cuda_model_pkg}" "${skala_model_urlpath}" "${skala_cuda_model_pkg}"
          fi
        fi
        echo "Installing skala CUDA model to ${pkg_install_dir}..."
        install -m 0644 "${skala_cuda_model_pkg}" "${pkg_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
      fi

      echo "Skala model installed successfully to ${pkg_install_dir}"
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage5/$(basename "${SCRIPT_NAME}")" \
          "${skala_model_pkg}" "${skala_cuda_model_pkg}"
      else
        write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage5/$(basename "${SCRIPT_NAME}")" \
          "${skala_model_pkg}"
      fi
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Skala from system paths ===================="
    check_dir "${SKALA_MODEL_ROOT}/share/skala/onedft_models"
    check_dir "${SKALA_MODEL_ROOT}/share/skala/onedft_models/${skala_model_pkg}"
    [ "${ENABLE_CUDA}" = "__TRUE__" ] && check_dir "${SKALA_MODEL_ROOT}/share/skala/onedft_models/${skala_cuda_model_pkg}"
    pkg_install_dir="${SKALA_MODEL_ROOT}"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Skala to custom path ===================="
    pkg_install_dir="${with_skala}"
    check_dir "${pkg_install_dir}/share/skala/onedft_models"
    check_dir "${pkg_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
    [ "${ENABLE_CUDA}" = "__TRUE__" ] && check_dir "${pkg_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
    ;;
esac

if [ "${with_skala}" != "__DONTUSE__" ]; then
  skala_model="${pkg_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
  if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
    skala_cuda_model="${pkg_install_dir}/share/skala/onedft_models/${skala_cuda_model_pkg}"
  fi
    cat << EOF > "${INSTALLDIR}/setup_skala"
export SKALA_VER="${skala_model_ver}"
export SKALA_ROOT="${pkg_install_dir}"
export SKALA_MODEL="${skala_model}"
EOF
   if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
     cat >> "${INSTALLDIR}/setup_skala" << EOF
export SKALA_CUDA_MODEL="${skala_cuda_model}"
EOF
   fi
   if [ "${with_skala}" != "__SYSTEM__" ]; then
     cat << EOF >> "${INSTALLDIR}/setup_skala"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
   fi
   filter_setup "${INSTALLDIR}/setup_skala" "${SETUPFILE}"
fi

load "${INSTALLDIR}/setup_skala"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "skala"
