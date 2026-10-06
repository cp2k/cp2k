#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

libtorch_ver="2.7.1"
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

    # Select a broadly available LibTorch/CUDA combination based on GPUVER:
    #   Pascal/Volta/Turing  2.5.1 + CUDA 11.8
    #   Ampere/Hopper        2.5.1 + CUDA 12.4
    #   Blackwell            2.7.1 + CUDA 12.8
    # CPU builds use the version mirrored on the CP2K download server.
    if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
      case "${GPUVER}" in
        P100 | V100)
          libtorch_ver="2.5.1"
          libtorch_cuda_suffix="cu118"
          libtorch_sha256="9b524e24c0ea15f191cfe6461594a3f1773d5c866b804bffd99882937a28081b"
          ;;
        A100 | A40 | H100)
          libtorch_ver="2.5.1"
          libtorch_cuda_suffix="cu124"
          libtorch_sha256="08552b1d13de1389d01c64aaa5e2c6c5092f922618b0d32ccff8ed9099d75735"
          ;;
        B200 | GB10)
          libtorch_ver="2.7.1"
          libtorch_cuda_suffix="cu128"
          libtorch_sha256="ae513b437ae99150744ef1d06b02a4ecbbb9275c9ffe540c88909623e3293041"
          ;;
        *)
          # Very old (K20X/K40/K80) or otherwise unlisted GPUs.
          libtorch_ver="2.5.1"
          libtorch_cuda_suffix="cu118"
          libtorch_sha256="9b524e24c0ea15f191cfe6461594a3f1773d5c866b804bffd99882937a28081b"
          ;;
      esac
      # Readable CUDA version, e.g. cu124 -> 12.4.
      libtorch_cuda_digits="${libtorch_cuda_suffix#cu}"
      libtorch_cuda_version="${libtorch_cuda_digits:0:2}.${libtorch_cuda_digits:2:1}"
      archive_file="libtorch-cxx11-abi-shared-with-deps-${libtorch_ver}+${libtorch_cuda_suffix}.zip"
      libtorch_url="https://download.pytorch.org/libtorch/${libtorch_cuda_suffix}/${archive_file}"
      echo "Installing CUDA-enabled libtorch ${libtorch_ver} (CUDA ${libtorch_cuda_version}) for ${GPUVER}"
    else
      libtorch_ver="2.7.1"
      libtorch_sha256="63d572598c8d532128a335018913e795c1bbb32602ce378896dc8cfbb5590976"
      archive_file="libtorch-cxx11-abi-shared-with-deps-${libtorch_ver}+cpu.zip"
      libtorch_url=""
      echo "Installing CPU-only libtorch ${libtorch_ver}"
    fi

    pkg_install_dir="${INSTALLDIR}/libtorch-${libtorch_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "libtorch-${libtorch_ver} is already installed, skipping it."
    else
      if [ -n "${libtorch_url}" ]; then
        echo "Downloading from: ${libtorch_url}"
        wget --quiet "${libtorch_url}" -O "${archive_file}" ||
          report_error ${LINENO} "Failed to download ${archive_file} from ${libtorch_url}"
      else
        retrieve_package "${libtorch_sha256}" "${archive_file}"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d libtorch ] && rm -rf libtorch
      [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
      unzip -q ${archive_file}
      mv libtorch ${pkg_install_dir}

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename "${SCRIPT_NAME}")"

      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        report_warning ${LINENO} \
          "A prebuilt CUDA-enabled libtorch (version ${libtorch_ver}, CUDA ${libtorch_cuda_version})
was installed. It is highly recommended to install your own libtorch matching your exact GPU
and CUDA version and point CP2K to it: --with-libtorch=system searches for it in the system
paths, or --with-libtorch=<path_to_libtorch> uses a libtorch you
built yourself. Otherwise poor performance is expected and the calculation may even fail."
      else
        report_warning ${LINENO} \
          "A prebuilt CPU-only libtorch (version ${libtorch_ver}) was installed. It is highly
recommended to install your own libtorch matching your platform and point CP2K to it:
--with-libtorch=system searches for it in the system paths, or
--with-libtorch=<path_to_libtorch> uses a libtorch you built yourself. Otherwise poor
performance is expected and the calculation may even fail."
      fi
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
