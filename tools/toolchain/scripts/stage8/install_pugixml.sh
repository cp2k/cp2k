#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

pugixml_ver="1.14"
pugixml_sha256="2f10e276870c64b1db6809050a75e11a897a8d7456c4be5c6b2e35a11168a015"
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
    pugixml_ROOT="${pkg_install_dir}"
    PUGIXML_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding pugixml from system paths ===================="
    check_command pkg-config --modversion pugixml
    add_lib_from_paths PUGIXML_LDFLAGS "libpugixml.*" $LIB_PATHS
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking pugixml to user paths ===================="
    pkg_install_dir="${with_pugixml}"

    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    PUGIXML_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && PUGIXML_LIBDIR="${pkg_install_dir}/lib64"

    check_dir "${PUGIXML_LIBDIR}"
    check_dir "${pkg_install_dir}/include/pugixml"
    PUGIXML_CFLAGS="-I'${pkg_install_dir}/include/pugixml'"
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
export PUGIXML_LIBS="-lpugixml"
export pugixml_ROOT="${pkg_install_dir}"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_pugixml"
export PUGIXML_LDFLAGS="${SPLA_LDFLAGS}"
export PUGIXML_LIBRARY="-lpugixml"
export pugixml_ROOT="$pkg_install_dir"
export PUGIXML_VERSION=${pugixml-ver}
export CP_LIBS="IF_MPI(${PUGIXML_LIBS}|) \${CP_LIBS}"
EOF
  cat << EOF >> "${BUILDDIR}/setup_pugixml"
export CP_LDFLAGS="\${CP_LDFLAGS} ${PUGIXML_LDFLAGS}"
EOF
  cat "${BUILDDIR}/setup_pugixml" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_pugixml"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "pugixml"
