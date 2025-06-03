#!/bin/bash -e

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

greenx_ver="2.2"
greenx_sha256="cf0abb77cc84a3381a690a6ac7ca839da0007bb9e6120f3f25e47de50e29431f"

# shellcheck disable=SC1091
source "${SCRIPT_DIR}"/common_vars.sh
# shellcheck disable=SC1091
source "${SCRIPT_DIR}"/tool_kit.sh
# shellcheck disable=SC1091
source "${SCRIPT_DIR}"/signal_trap.sh
# shellcheck disable=SC1091
source "${INSTALLDIR}"/toolchain.conf
# shellcheck disable=SC1091
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_greenx" ] && rm "${BUILDDIR}/setup_greenx"

GREENX_CFLAGS=""
GREENX_LDFLAGS=""
GREENX_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
with_greenx=${with_greenx:__DONTUSE__}
with_gmp=${with_greenx:__DONTUSE__}
case "$with_greenx" in
  __INSTALL__)
    echo "==================== Installing GreenX ===================="
    pkg_install_dir="${INSTALLDIR}/greenX-${greenx_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "greenX-${greenx_ver} is already installed, skipping it."
    else
      if [ -f greenX-${greenx_ver}.tar.gz ]; then
        echo "greenX-${greenx_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${greenx_sha256}" "greenX-${greenx_ver}.tar.gz"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d greenX-${greenx_ver} ] && rm -rf greenX-${greenx_ver}
      tar -xzf greenX-${greenx_ver}.tar.gz
      cd greenX-${greenx_ver}
      mkdir build
      cd build
      # CP2K only uses the AC and Minimax components
      # TODO : Add the stdc++ to libs when libint not used?
      if [ "$with_gmp" != "__DONTUSE__" ]; then
        gmp_flag=ON
      else
        gmp_flag=OFF
      fi
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DAC_COMPONENT=ON \
        -DMINIMAX_COMPONENT=ON \
        -DENABLE_GNU_GMP=$gmp_flag \
        -DENABLE_GREENX_CTEST=OFF \
        .. > configure.log 2>&1 || tail -n "${LOG_LINES}" configure.log
      make -j "$(get_nprocs)" > make.log 2>&1 || tail -n "${LOG_LINES}" make.log
      make install > install.log 2>&1 || tail -n "${LOG_LINES}" install.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage3/$(basename "${SCRIPT_NAME}")"
      cd ..
    fi
    GREENX_CFLAGS="-I'${pkg_install_dir}/include/modules'"
    GREENX_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding GreenX from system paths ===================="
    check_lib -lGXCommon "greenx"
    check_lib -lgx_ac "greenx"
    check_lib -lgx_minimax "greenx"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking GreenX to user paths ===================="
    pkg_install_dir="$with_greenx"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include/modules"
    GREENX_CFLAGS="-I'${pkg_install_dir}/include/modules'"
    GREENX_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_greenx" != "__DONTUSE__" ]; then
  GREENX_LIBS=" -lGXCommon -lgx_ac -lgx_minimax"
  cat << EOF > "${BUILDDIR}/setup_greenx"
export GREENX_VER="${greenx_ver}"
EOF
  if [ "$with_greenx" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_greenx"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_greenx"
export GREENX_CFLAGS="${GREENX_CFLAGS}"
export GREENX_LDFLAGS="${GREENX_LDFLAGS}"
export GREENX_LIBS="${GREENX_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__GREENX"
export CP_CFLAGS="\${CP_CFLAGS} ${GREENX_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${GREENX_LDFLAGS}"
export CP_LIBS="${GREENX_LIBS} \${CP_LIBS}"
export GREENX_ROOT="${pkg_install_dir}"
EOF
  cat "${BUILDDIR}/setup_greenx" >> "$SETUPFILE"
fi

load "${BUILDDIR}/setup_greenx"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "greenx"
