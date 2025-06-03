#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

dftd4_ver="3.7.0"
dftd4_sha256="2e0d3504038358b8a82fdd21912b7765d416a58ebedbdd44f2ca8d2e88339ad7"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_dftd4" ] && rm "${BUILDDIR}/setup_dftd4"

DFTD4_DFLAGS=''
DFTD4_CFLAGS=''
DFTD4_LDFLAGS=''
DFTD4_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_dftd4" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing DFTD4 ===================="
    require_env OPENBLAS_ROOT
    require_env MATH_LIBS

    pkg_install_dir="${INSTALLDIR}/dftd4-${dftd4_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "dftd4_dist-${dftd4_ver} is already installed, skipping it."
    else
      if [ -f dftd4-${dftd4_ver}.tar.gz ]; then
        echo " dftd4-${dftd4_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${dftd4_sha256}" "dftd4-${dftd4_ver}.tar.gz"
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d dftd4-${dftd4_ver} ] && rm -rf dftd4-${dftd4_ver}
      tar -xzf dftd4-${dftd4_ver}.tar.gz
      cd dftd4-${dftd4_ver}/subprojects/dftd4

      rm -Rf build
      mkdir build
      cd build

      CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${OPENBLAS_ROOT}" cmake \
        -B . -G Ninja \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      cmake --build . -j $(get_nprocs) >> build.log 2>&1 || tail -n ${LOG_LINES} build.log
      cmake --install . >> install.log 2>&1 || tail -n ${LOG_LINES} install.log

      cd ..
    fi
    write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    ;;

  __SYSTEM__)
    echo "==================== Finding DFTD4 from system paths ===================="
    check_command pkg-config --modversion dftd4
    add_include_from_paths DFTD4_CFLAGS "dftd4.h" $INCLUDE_PATHS
    add_include_from_paths DFTD4_CFLAGS "dftd4.mod" $INCLUDE_PATHS
    add_include_from_paths DFTD4_CFLAGS "mctc_io.mod" $INCLUDE_PATHS
    add_lib_from_paths DFTD4_LDFLAGS "libdftd4.*" $LIB_PATHS
    ;;

  *)
    echo "==================== Linking DFTD4 to user paths ===================="
    pkg_install_dir="$with_dftd4"
    check_dir "${pkg_install_dir}/include"
    ;;

esac

if [ "$with_dftd4" != "__DONTUSE__" ]; then

  DFTD4_DFLAGS="-D__DFTD4"
  DFTD4_LIBS="-ldftd4 -lmstore -lmulticharge -lmctc-lib"

  cat << EOF > "${BUILDDIR}/setup_dftd4"
export DFTD4_VER="${dftd4_ver}"
EOF

  if [ "$with_dftd4" != "__SYSTEM__" ]; then
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "multicharge.mod")
    MCHARGE=${DFTD4_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "mstore.mod")
    MSTORE=${DFTD4_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "mctc_io.mod")
    MCTC=${DFTD4_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "dftd4.mod")
    DFTD4=${DFTD4_LOC%/*}

    DFTD4_INCLUDE_DIRS="$pkg_install_dir/include"
    DFTD4_LINK_LIBRARIES="${pkg_install_dir}/lib"

    DFTD4_CFLAGS="-I'${MCHARGE}' -I'${MCTC}' -I'${DFTD4}'"
    DFTD4_LDFLAGS="-L'${DFTD4_LINK_LIBRARIES}' -Wl,-rpath,'${DFTD4_LINK_LIBRARIES}'"

    cat << EOF >> "${BUILDDIR}/setup_dftd4"
prepend_path LD_LIBRARY_PATH "${DFTD4_LINK_LIBRARIES}"
prepend_path LD_RUN_PATH "${DFTD4_LINK_LIBRARIES}"
prepend_path LIBRARY_PATH "${DFTD4_LINK_LIBRARIES}"
prepend_path CPATH "${DFTD4_INCLUDE_DIRS}"
prepend_path PKG_CONFIG_PATH "${DFTD4_LINK_LIBRARIES}/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi

  cat << EOF >> "${BUILDDIR}/setup_dftd4"
export MCHARGE="${MCHARGE}"
export MSTORE="${MSTORE}"
export MCTC="${MCTC}"
export DFTD4="${DFTD4}"
export DFTD4_INCLUDE_DIRS="${DFTD4_INCLUDE_DIRS}"
export DFTD4_LINK_LIBRARIES="${DFTD4_LINK_LIBRARIES}"
export DFTD4_ROOT="${pkg_install_dir}"
export DFTD4_DFLAGS="${DFTD4_DFLAGS}" 
export DFTD4_CFLAGS="${DFTD4_CFLAGS}"
export DFTD4_LDFLAGS="${DFTD4_LDFLAGS}"
export DFTD4_LIBS="${DFTD4_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} \${DFTD4_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} \${DFTD4_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} \${DFTD4_LDFLAGS}"
export CP_LIBS="\${DFTD4_LIBS} \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_dftd4" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_dftd4"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "dftd4"
