#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

gauxc_ver="1.0.0-skala-cp2k-rks-density-fix"
gauxc_rev="f6ab248a32108036c6e92886d78736d99fcdfc87"
gauxc_pkg="GauXC-${gauxc_rev}.tar.gz"
gauxc_sha256="6b2aca6fe7a498665f8d8e6f7f76c0ef1051d8c35bb23b58468226f7171a0376"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_gauxc" ] && rm "${BUILDDIR}/setup_gauxc"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

retrieve_github_archive() {
  local __sha256="$1"
  local __filename="$2"
  local __urlpath="$3"
  local __outfile="$4"
  if ! [ -f "${__outfile}" ]; then
    download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}"
  elif ! checksum "${__sha256}" "${__outfile}"; then
    echo "${__outfile} is found but checksum is wrong; delete and re-download"
    rm -vf "${__outfile}"
    download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}"
  else
    echo "${__outfile} is found and checksum is right"
  fi
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
      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf "GauXC-${gauxc_rev}" "${pkg_install_dir}"
      tar -xzf "${gauxc_pkg}"
      cd "GauXC-${gauxc_rev}"
      patch -l -p1 < "${SCRIPT_DIR}/stage3/gauxc-${gauxc_ver}.patch" \
        > gauxc_cp2k_rks_density_fix.patch.log 2>&1 || tail_excerpt gauxc_cp2k_rks_density_fix.patch.log
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

      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=OFF \
        -DGAUXC_ENABLE_C=ON \
        -DGAUXC_ENABLE_FORTRAN=ON \
        -DGAUXC_ENABLE_HOST=ON \
        -DGAUXC_ENABLE_MPI="${gauxc_enable_mpi}" \
        -DGAUXC_ENABLE_OPENMP="${gauxc_enable_openmp}" \
        -DGAUXC_ENABLE_TESTS=OFF \
        -DGAUXC_ENABLE_CUDA=OFF \
        -DGAUXC_ENABLE_HIP=OFF \
        -DGAUXC_ENABLE_HDF5=OFF \
        -DGAUXC_ENABLE_ONEDFT=ON \
        -DGAUXC_ENABLE_MAGMA=OFF \
        -DGAUXC_ENABLE_NCCL=OFF \
        -DGAUXC_ENABLE_CUTLASS=OFF \
        .. > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename ${SCRIPT_NAME})" \
        "${SCRIPT_DIR}/stage3/gauxc-${gauxc_ver}.patch" "${BUILDDIR}/${gauxc_pkg}"
    fi
    GAUXC_CFLAGS="-I'${pkg_install_dir}/include' -I'${pkg_install_dir}/include/gauxc/modules'"
    GAUXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding GauXC from system paths ===================="
    check_lib -lgauxc "GauXC"
    add_include_from_paths -p GAUXC_CFLAGS "gauxc" $INCLUDE_PATHS
    add_include_from_paths GAUXC_CFLAGS "gauxc/modules/gauxc_status.mod" $INCLUDE_PATHS
    add_lib_from_paths GAUXC_LDFLAGS "libgauxc.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking GauXC to user paths ===================="
    pkg_install_dir="${with_gauxc}"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/include/gauxc/modules"
    GAUXC_CFLAGS="-I'${pkg_install_dir}/include' -I'${pkg_install_dir}/include/gauxc/modules'"
    GAUXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac

if [ "${with_gauxc}" != "__DONTUSE__" ]; then
  GAUXC_LIBS="-lgauxc -lintegratorxx -lexchcxx"
  GAUXC_DFLAGS="-D__GAUXC"
  gauxc_config_file="$(find_in_paths "gauxc/gauxc_config.f" $INCLUDE_PATHS)"
  if { [ "${gauxc_config_file}" != "__FALSE__" ] && grep -q "GAUXC_HAS_MPI" "${gauxc_config_file}"; } ||
    { [ -f "${pkg_install_dir}/include/gauxc/gauxc_config.f" ] &&
      grep -q "GAUXC_HAS_MPI" "${pkg_install_dir}/include/gauxc/gauxc_config.f"; }; then
    GAUXC_DFLAGS="${GAUXC_DFLAGS} -DGAUXC_HAS_MPI"
  fi

  cat << EOF > "${BUILDDIR}/setup_gauxc"
export GAUXC_VER="${gauxc_ver}"
export GAUXC_CFLAGS="${GAUXC_CFLAGS}"
export GAUXC_LDFLAGS="${GAUXC_LDFLAGS}"
export GAUXC_LIBS="${GAUXC_LIBS}"
export GAUXC_ROOT="${pkg_install_dir}"
EOF
  if [ "${with_gauxc}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_gauxc"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include/gauxc/modules"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_gauxc"
export CP_DFLAGS="\${CP_DFLAGS} ${GAUXC_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${GAUXC_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${GAUXC_LDFLAGS}"
export CP_LIBS="${GAUXC_LIBS} \${CP_LIBS}"
EOF
  filter_setup "${BUILDDIR}/setup_gauxc" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_gauxc"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "gauxc"
