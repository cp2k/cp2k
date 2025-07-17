#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

pugixml_ver="1.15"
pugixml_sha256="655ade57fa703fb421c2eb9a0113b5064bddb145d415dd1f88c79353d90d511a"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_pugixml" ] && rm "${BUILDDIR}/setup_pugixml"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_pugixml}" in
  __INSTALL__)
    echo "==================== Installing pugixml ===================="
    pkg_install_dir="${INSTALLDIR}/pugixml-${pugixml_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "pugixml-${pugixml_ver} is already installed, skipping it."
    else
      if [ -f pugixml-${pugixml_ver}.tar.gz ]; then
        echo "pugixml-${pugixml_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${pugixml_sha256}" "pugixml-${pugixml_ver}.tar.gz"

      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d pugixml-${pugixml_ver} ] && rm -rf pugixml-${pugixml_ver}
      tar -xzf pugixml-${pugixml_ver}.tar.gz
      cd pugixml-${pugixml_ver}
      mkdir -p build
      cd build
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    PUGIXML_CFLAGS="-I'${pkg_install_dir}/include'"
    PUGIXML_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding pugixml from system paths ===================="
    check_command pkg-config --modversion pugixml
    add_include_from_paths PUGIXML_CFLAGS "pugixml.hpp" ${INCLUDE_PATHS}
    add_lib_from_paths PUGIXML_LDFLAGS "libpugixml.*" ${LIB_PATHS}
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking pugixml to user paths ===================="
    pkg_install_dir="${with_pugixml}"
    # Use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    PUGIXML_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && PUGIXML_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${PUGIXML_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    PUGIXML_CFLAGS="-I'${pkg_install_dir}/include'"
    PUGIXML_LDFLAGS="-L'${PUGIXML_LIBDIR}' -Wl,-rpath,'${PUGIXML_LIBDIR}'"
    ;;
esac

if [ "${with_pugixml}" != "__DONTUSE__" ]; then
  PUGIXML_LIBS="-lpugixml"
  cat << EOF > "${BUILDDIR}/setup_pugixml"
export PUGIXML_VER="${pugixml_ver}"
EOF
  if [ "${with_pugixml}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_pugixml"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CPATH "${pkg_install_dir}/include"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_pugixml"
export PUGIXML_CFLAGS="${PUGIXML_CFLAGS}"
export PUGIXML_LDFLAGS="${PUGIXML_LDFLAGS}"
export PUGIXML_LIBRARY="-lpugixml"
export CP_CFLAGS="\${CP_CFLAGS} ${PUGIXML_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${PUGIXML_LDFLAGS}"
export CP_LIBS="IF_MPI(${PUGIXML_LIBS}|) \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_pugixml" >> ${SETUPFILE}
fi

load "${BUILDDIR}/setup_pugixml"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "pugixml"
