#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

gauxc_ver="1.1-skala-cp2k-fixes"
gauxc_rev="ea0f1fc7c02497aa950376eda61686ba30256999"
gauxc_pkg="GauXC-${gauxc_rev}.tar.gz"
gauxc_sha256="10b89536abd8d43b0f10ae33b53c8218c45f23a375bbd9783b62e8e6e606b579"
nlohmann_json_pkg="nlohmann-json-3.12.0-include.zip"
nlohmann_json_sha256="b8cb0ef2dd7f57f18933997c9934bb1fa962594f701cd5a8d3c2c80541559372"
nlohmann_json_urlpath="https://github.com/nlohmann/json/releases/download/v3.12.0"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_gauxc" ] && rm "${BUILDDIR}/setup_gauxc"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

source "${INSTALLDIR}/setup_skala"

retrieve_github_archive() {
  local __sha256="$1" __filename="$2" __urlpath="$3" __outfile="$4" __attempt
  if ! [ -f "${__outfile}" ] || ! checksum "${__sha256}" "${__outfile}"; then
    if [ -f "${__outfile}" ]; then
      echo "${__outfile} checksum wrong, deleting..."
      rm -vf "${__outfile}"
    fi
    for __attempt in 1 2 3 4 5; do
      if download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}"; then
        return
      fi
      echo "Download attempt ${__attempt} for ${__filename} failed."
      sleep $((__attempt * 5))
    done
    return 1
  fi
  echo "${__outfile} is found and checksum is right"
}

case "${with_gauxc}" in
  __INSTALL__)
    echo "==================== Installing GauXC ===================="
    pkg_install_dir="${INSTALLDIR}/gauxc-${gauxc_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "gauxc-${gauxc_ver} is already installed, skipping it."
    else
      retrieve_github_archive "${gauxc_sha256}" "${gauxc_rev}.tar.gz" \
        "https://github.com/wavefunction91/GauXC/archive" \
        "${gauxc_pkg}"
      retrieve_github_archive "${nlohmann_json_sha256}" "include.zip" \
        "${nlohmann_json_urlpath}" \
        "${nlohmann_json_pkg}"
      # Skala model is installed by shared install_skala.sh script
      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf "GauXC-${gauxc_rev}" "${pkg_install_dir}"
      tar -xzf "${gauxc_pkg}"
      cd "GauXC-${gauxc_rev}"

      # CMake configuration modification
      if ! patch -l -p1 < "${SCRIPT_DIR}/stage6/gauxc-${gauxc_ver}.patch" \
        > gauxc_cp2k_rks_density_fix.patch.log 2>&1; then
        tail_excerpt gauxc_cp2k_rks_density_fix.patch.log
      fi
      install -m 0644 "${SCRIPT_DIR}/stage6/exchcxx-disable-builtin.patch" \
        cmake/exchcxx-disable-builtin.patch
      if ! patch -l -p1 < "${SCRIPT_DIR}/stage6/gauxc-libxc-only-exchcxx.patch" \
        > gauxc_libxc_only_exchcxx.patch.log 2>&1; then
        tail_excerpt gauxc_libxc_only_exchcxx.patch.log
      fi

      sed -i.bak '/find_dependency.*nlohmann_json/d' cmake/gauxc-config.cmake.in
      rm -f cmake/gauxc-config.cmake.in.bak

      # Regenerate Obara-Saika kernels with lower unroll threshold
      cd src/xc_integrator/local_work_driver/host/obara_saika/generator
      make > make.log 2>&1 || tail_excerpt make.log
      cd ../src
      ../generator/generate_cpu_code.x 4 4
      perl -pi -e 's#\.\./include/#../include/cpu/#g' \
        integral_*.cxx \
        integral_*.hpp \
        obara_saika_integrals.cxx
      cd "${BUILDDIR}/GauXC-${gauxc_rev}"

      mkdir -p build
      cd build

      if [ "${MPI_MODE}" = "no" ]; then
        gauxc_enable_mpi="OFF"
      else
        gauxc_enable_mpi="ON"
      fi
      if [ "${ENABLE_OMP}" = "__TRUE__" ]; then
        gauxc_enable_openmp="ON"
      else
        gauxc_enable_openmp="OFF"
      fi
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        gauxc_enable_cuda="ON"
        gauxc_cuda_architectures_option="-DCMAKE_CUDA_ARCHITECTURES=${ARCH_NUM}"
      else
        gauxc_enable_cuda="OFF"
        gauxc_cuda_architectures_option=""
      fi
      if [ "${ENABLE_GAUXC_CUTLASS}" = "__TRUE__" ]; then
        gauxc_enable_cutlass="ON"
      else
        gauxc_enable_cutlass="OFF"
      fi

      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DGAUXC_ENABLE_C=ON \
        -DGAUXC_ENABLE_FORTRAN=ON \
        -DGAUXC_ENABLE_HOST=ON \
        -DGAUXC_ENABLE_MPI="${gauxc_enable_mpi}" \
        -DGAUXC_ENABLE_OPENMP="${gauxc_enable_openmp}" \
        -DGAUXC_ENABLE_TESTS=OFF \
        -DGAUXC_ENABLE_CUDA="${gauxc_enable_cuda}" \
        -DGAUXC_ENABLE_HIP=OFF \
        -DGAUXC_ENABLE_HDF5=OFF \
        -DGAUXC_ENABLE_ONEDFT=ON \
        -DGAUXC_ENABLE_EXCHCXX_BUILTIN=OFF \
        -DGAUXC_NLOHMANN_JSON_URL="${BUILDDIR}/${nlohmann_json_pkg}" \
        -DGAUXC_ENABLE_MAGMA=OFF \
        -DGAUXC_ENABLE_NCCL=OFF \
        -DGAUXC_ENABLE_CUTLASS="${gauxc_enable_cutlass}" \
        ${gauxc_cuda_architectures_option} \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename "${SCRIPT_NAME}")" \
        "${SCRIPT_DIR}/stage6/gauxc-${gauxc_ver}.patch" \
        "${SCRIPT_DIR}/stage6/gauxc-libxc-only-exchcxx.patch" \
        "${SCRIPT_DIR}/stage6/exchcxx-disable-builtin.patch" "${BUILDDIR}/${gauxc_pkg}" \
        "${BUILDDIR}/${nlohmann_json_pkg}"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding GauXC from system paths ===================="
    check_lib -lgauxc "GauXC"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libgauxc.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking GauXC to user paths ===================="
    pkg_install_dir="${with_gauxc}"
    GAUXC_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && GAUXC_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${GAUXC_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/include/gauxc/modules"
    ;;
esac

if [ "${with_gauxc}" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_gauxc"
export GAUXC_VER="${gauxc_ver}"
export GAUXC_ROOT="${pkg_install_dir}"
export SKALA_MODEL="${SKALA_MODEL:-}"
export SKALA_CUDA_MODEL="${SKALA_CUDA_MODEL:-}"
EOF
  if [ "${with_gauxc}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_gauxc"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_gauxc" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_gauxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "gauxc"
