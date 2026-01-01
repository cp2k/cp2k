#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

mcl_ver="3.0.0"
mcl_sha256="3e740582836fe90e04a693cfc5a219826bcac03217f70ea5570bad6aeafda685"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_mcl" ] && rm "${BUILDDIR}/setup_mcl"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_mcl:=__INSTALL__}" in
  __INSTALL__)
    echo "==================== Installing MCL ===================="
    pkg_install_dir="${INSTALLDIR}/mcl-${mcl_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libmcl-${mcl_ver} is already installed, skipping it."
    else
      if [ -f mcl-${mcl_ver}.tar.gz ]; then
        echo "mcl-${mcl_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${mcl_sha256}" mcl-${mcl_ver}.tar.gz
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d mcl-${mcl_ver} ] && rm -rf mcl-${mcl_ver}
      tar -xzf mcl-${mcl_ver}.tar.gz

      mkdir "mcl-${mcl_ver}/build"
      cd "mcl-${mcl_ver}/build"

      cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=YES \
        -DBUILD_FORTRAN_API=YES \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . > build.log 2>&1 || tail -n ${LOG_LINES} build.log
      CMAKE_BUILD_PARALLEL_LEVEL="$(get_nprocs)" cmake --build . --target install > install.log 2>&1 || tail -n ${LOG_LINES} install.log

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename "${SCRIPT_NAME}")"
    fi
    MCL_LDFLAGS="-L'${pkg_install_dir}/lib/MiMiC' -Wl,-rpath,'${pkg_install_dir}/lib/MiMiC'"
    ;;
  __SYSTEM__)
    echo "==================== Finding MCL from system paths ===================="
    check_lib -lmclf "mcl"
    add_lib_from_paths MCL_LDFLAGS "mclf.*" "$LIB_PATHS"
    add_lib_from_paths MCL_LDFLAGS "mcl.*" "$LIB_PATHS"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking MCL to user paths ===================="
    pkg_install_dir="${with_mcl}"

    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    MCL_LIBDIR="${pkg_install_dir}/lib/MiMiC"
    [ -d "${pkg_install_dir}/lib64" ] && MCL_LIBDIR="${pkg_install_dir}/lib64/MiMiC"

    check_dir "${MCL_LIBDIR}"
    MCL_LDFLAGS="-L'${MCL_LIBDIR}' -Wl,-rpath,'$MCL_LIBDIR}'"
    ;;
esac

if [ "$with_mcl" != "__DONTUSE__" ]; then
  MCL_LIBS="-lmcl -lmclf -lstdc++"
  if [ "$with_mcl" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_mcl"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib/MiMiC"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib/MiMiC"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib/MiMiC"
export MCL_LIBS="${MCL_LIBS}"
export MCL_ROOT="${pkg_install_dir}"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_mcl"
export MCL_VER="${mcl_ver}"
export MCL_ROOT="${pkg_install_dir}"
export MCL_LDFLAGS="${MCL_LDFLAGS}"
export MCL_LIBS="${MCL_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__MIMIC"
export CP_LDFLAGS="\${CP_LDFLAGS} ${MCL_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${MCL_LIBS}"
EOF
  cat "${BUILDDIR}/setup_mcl" >> "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_mcl"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "MCL"
