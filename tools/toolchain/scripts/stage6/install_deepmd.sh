#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

deepmd_ver="2.2.6"
deepmd_sha256="fd9aa49194e0030b312c4741a3bd26b94d06110abc76ffdef7a6d236286922c2"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_deepmd" ] && rm "${BUILDDIR}/setup_deepmd"

DEEPMD_LDFLAGS=''
DEEPMD_LIBS=''
DEEPMD_CXXFLAGS='-std=gnu++11 '

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_deepmd" in
  __INSTALL__)
    echo "==================== Installing DeePMD ===================="
    pkg_install_dir="${INSTALLDIR}/libdeepmd_c-${deepmd_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    deepmd_root="${pkg_install_dir}"
    if verify_checksums "${install_lock_file}"; then
      echo "libdeepmd_c-${deepmd_ver} is already installed, skipping it."
    else
      if [ -f deepmd-kit-${deepmd_ver}.tar.gz ]; then
        echo "libdeepmd_c-${deepmd_ver}.tar.gz is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} -o libdeepmd_c-${deepmd_ver}.tar.gz ${deepmd_sha256} \
          https://github.com/deepmodeling/deepmd-kit/releases/download/v${deepmd_ver}/libdeepmd_c.tar.gz
      fi
      [ -d libdeepmd_c ] && rm -rf libdeepmd_c
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf libdeepmd_c-${deepmd_ver}.tar.gz
      mv libdeepmd_c ${pkg_install_dir}
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    DEEPMD_DFLAGS="-D__DEEPMD"
    DEEPMD_CFLAGS="-I'${pkg_install_dir}/include'"
    if [ "$DEEPMD_MODE" == "cpu" ]; then
      DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_op -ldeepmd_c -Wl,-rpath='${pkg_install_dir}/lib'"
    elif [ "$DEEPMD_MODE" == "cuda" ]; then
      DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_op -ldeepmd_c -ldeepmd_dyn_cudart -ldeepmd_op_cuda -Wl,-rpath='${pkg_install_dir}/lib'"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding DeePMD from system paths ===================="
    check_lib -ldeepmd "DEEPMD"
    add_lib_from_paths DEEPMD_LDFLAGS "libdeepmd*" $LIB_PATHS
    add_include_from_paths DEEPMD_CFLAGS "deepmd" $INCLUDE_PATHS
    add_include_from_paths DEEPMD_CXXFLAGS "deepmd" $INCLUDE_PATHS
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
    if [ "$DEEPMD_MODE" == "cpu" ]; then
      DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_op -ldeepmd_c -Wl,-rpath='${pkg_install_dir}/lib'"
    elif [ "$DEEPMD_MODE" == "cuda" ]; then
      DEEPMD_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,--no-as-needed -ldeepmd_op -ldeepmd_c -ldeepmd_dyn_cudart -ldeepmd_op_cuda -Wl,-rpath='${pkg_install_dir}/lib'"
    fi
  ;;
esac

if [ "$with_deepmd" != "__DONTUSE__" ]; then
  if [ "$DEEPMD_MODE" == "cpu" ]; then
    DEEPMD_LIBS='-ldeepmd_op -ldeepmd_op -ldeepmd_c -lstdc++'
  elif [ "$DEEPMD_MODE" == "cuda" ]; then
    DEEPMD_LIBS='-ldeepmd_op -ldeepmd_op -ldeepmd_c -ldeepmd_dyn_cudart -ldeepmd_op_cuda -lstdc++'
  fi
  if [ "$with_deepmd" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_deepmd"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
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
fi

load "${BUILDDIR}/setup_deepmd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "deepmd"
