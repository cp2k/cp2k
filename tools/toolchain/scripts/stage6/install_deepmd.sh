#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

deepmd_ver="3.0.1"
deepmd_pkg="deepmd-kit-${deepmd_ver}.tar.gz"
deepmd_sha256="e842edbc2714bc948ce708c411e5fed751e67c88d5c493c2978f11c849027dca"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_deepmd" ] && rm "${BUILDDIR}/setup_deepmd"

DEEPMD_LDFLAGS=''
DEEPMD_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_deepmd" in
  __INSTALL__)
    echo "==================== Installing DeePMD ===================="
    pkg_install_dir="${INSTALLDIR}/deepmd-kit-${deepmd_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    deepmd_root="${pkg_install_dir}"
    if verify_checksums "${install_lock_file}"; then
      echo "libdeepmd_c-${deepmd_ver} is already installed, skipping it."
    else
      if [ -f ${deepmd_pkg} ]; then
        echo "${deepmd_pkg} is found"
      else
        download_pkg_from_cp2k_org "${deepmd_sha256}" "${deepmd_pkg}"
      fi
      [ -d deepmd-kit-${deepmd_ver} ] && rm -rf deepmd-kit-${deepmd_ver}
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf ${deepmd_pkg}
      cd deepmd-kit-${deepmd_ver}/source
      # Workaround for https://github.com/deepmodeling/deepmd-kit/issues/4569
      sed -i /CXX_STANDARD/d CMakeLists.txt

      # PR 4577: https://github.com/deepmodeling/deepmd-kit/pull/4577
      patch -p2 CMakeLists.txt < ${SCRIPT_DIR}/stage6/deepmd-kit_4577.patch

      mkdir build
      cd build
      cmake \
        -DENABLE_PYTORCH=TRUE \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_CXX_STANDARD=17 \
        -DCMAKE_CXX_STANDARD_REQUIRED=TRUE \
        .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j deepmd_c > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    DEEPMD_DFLAGS="-D__DEEPMD"
    DEEPMD_CFLAGS="-I'${pkg_install_dir}/include'"
    DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_c -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding DeePMD from system paths ===================="
    check_lib -ldeepmd "DEEPMD"
    add_lib_from_paths DEEPMD_LDFLAGS "libdeepmd*" $LIB_PATHS
    add_include_from_paths DEEPMD_CFLAGS "deepmd" $INCLUDE_PATHS
    DEEPMD_DFLAGS="-D__DEEPMD"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking DEEPMD to user paths ===================="
    pkg_install_dir="$with_deepmd"
    check_dir "${pkg_install_dir}/include/deepmd"
    check_dir "${pkg_install_dir}/lib"
    DEEPMD_DFLAGS="-D__DEEPMD"
    DEEPMD_CFLAGS="-I'${pkg_install_dir}/include'"
    DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_c -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_deepmd" != "__DONTUSE__" ]; then
  DEEPMD_LIBS='-ldeepmd_c -lstdc++'
  cat << EOF > "${BUILDDIR}/setup_deepmd"
export DEEPMD_VER="${deepmd_ver}"
EOF
  if [ "$with_deepmd" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_deepmd"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
    cat "${BUILDDIR}/setup_deepmd" >> $SETUPFILE
  fi

  cat << EOF >> "${BUILDDIR}/setup_deepmd"
export DEEPMD_DFLAGS="${DEEPMD_DFLAGS}"
export DEEPMD_CFLAGS="${DEEPMD_CFLAGS}"
export DEEPMD_LDFLAGS="${DEEPMD_LDFLAGS}"
export DEEPMD_LIBS="${DEEPMD_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} ${DEEPMD_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${DEEPMD_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${DEEPMD_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${DEEPMD_LIBS}"
EOF
  cat << EOF >> "${INSTALLDIR}/lsan.supp"
# leaks related to DeePMD-kit and TensorFlow
leak:DP_NewDeepPot
leak:deepmd::AtomMap::AtomMap
EOF
fi

load "${BUILDDIR}/setup_deepmd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "deepmd"
