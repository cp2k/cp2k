#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

ftorch_ver="1.0.0"
ftorch_sha256="c4b6741e582623b7ecaecd59d02f779e8a6f6017f8068c85da8a034f468df375"
ftorch_pkg="FTorch-${ftorch_ver}.tar.gz"
ftorch_urlpath="https://github.com/Cambridge-ICCS/FTorch/archive/refs/tags"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_ftorch" ] && rm "${BUILDDIR}/setup_ftorch"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

source "${BUILDDIR}/setup_skala"

retrieve_github_archive() {
  local __sha256="$1"
  local __filename="$2"
  local __urlpath="$3"
  local __outfile="$4"
  local __attempt
  if ! [ -f "${__outfile}" ]; then
    for __attempt in 1 2 3 4 5; do
      download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}" && return
      echo "Download attempt ${__attempt} for ${__filename} failed."
      sleep $((__attempt * 5))
    done
    return 1
  elif ! checksum "${__sha256}" "${__outfile}"; then
    echo "${__outfile} is found but checksum is wrong; delete and re-download"
    rm -vf "${__outfile}"
    for __attempt in 1 2 3 4 5; do
      download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}" && return
      echo "Download attempt ${__attempt} for ${__filename} failed."
      sleep $((__attempt * 5))
    done
    return 1
  else
    echo "${__outfile} is found and checksum is right"
  fi
}

case "${with_skala_ftorch}" in
  __INSTALL__)
    echo "==================== Installing FTorch + skala bindings ===================="
    pkg_install_dir="${INSTALLDIR}/ftorch-${ftorch_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "ftorch-${ftorch_ver} is already installed, skipping it."
    else
      retrieve_github_archive "${ftorch_sha256}" "v${ftorch_ver}.tar.gz" \
        "${ftorch_urlpath}" "${ftorch_pkg}"
      # Skala model is now installed by shared install_skala.sh script

      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf "FTorch-${ftorch_ver}" "${pkg_install_dir}"
      tar -xzf "${ftorch_pkg}"

      mkdir -p build_ftorch
      cd build_ftorch

      cmake "${BUILDDIR}/FTorch-${ftorch_ver}" \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=ON \
        > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log

      cd "${BUILDDIR}"

      write_checksums "${install_lock_file}" \
        "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})" \
        "${BUILDDIR}/${ftorch_pkg}"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding FTorch from system paths ===================="
    check_lib -lftorch "FTorch"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libftorch.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking skala-ftorch to user paths ===================="
    pkg_install_dir="${with_skala_ftorch}"
    FTORCH_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && FTORCH_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${FTORCH_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac

if [ "${with_skala_ftorch}" != "__DONTUSE__" ]; then
  if [ -n "${pkg_install_dir:-}" ]; then
    ftorch_skala_model="${SKALA_MODEL}"
  else
    ftorch_skala_model=""
  fi
  cat << EOF > "${BUILDDIR}/setup_ftorch"
export FTORCH_VER="${ftorch_ver}"
export FTORCH_ROOT="${pkg_install_dir}"
export SKALA_MODEL="${ftorch_skala_model}"
EOF
  if [ "${with_skala_ftorch}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_ftorch"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_ftorch" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_ftorch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "ftorch"
