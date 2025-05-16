#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

tblite_ver="0.4.0"
tblite_sha256="c4a67dfbe04827095fd7598183e076fa3017a5a475c4f90fd28e78992dc19ea7"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_tblite" ] && rm "${BUILDDIR}/setup_tblite"

TBLITE_DFLAGS=''
TBLITE_CFLAGS=''
TBLITE_LDFLAGS=''
TBLITE_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_tblite" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing tblite ===================="
    require_env OPENBLAS_ROOT
    require_env MATH_LIBS

    pkg_install_dir="${INSTALLDIR}/tblite-${tblite_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "tblite-${tblite_ver} is already installed, skipping it."
    else
      if [ -f tblite-${tblite_ver}.tar.gz ]; then
        echo "tblite-${tblite_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${tblite_sha256}" "tblite-${tblite_ver}.tar.gz"
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d tblite-${tblite_ver} ] && rm -rf tblite-${tblite_ver}
      tar -xzf tblite-${tblite_ver}.tar.gz
      cd tblite-${tblite_ver}

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
    echo "==================== Finding tblite from system paths ===================="
    check_command pkg-config --modversion tblite
    add_include_from_paths TBLITE_CFLAGS "tblite.h" $INCLUDE_PATHS
    add_include_from_paths TBLITE_CFLAGS "tblite.mod" $INCLUDE_PATHS
    add_include_from_paths TBLITE_CFLAGS "mctc_io.mod" $INCLUDE_PATHS
    add_lib_from_paths TBLITE_LDFLAGS "libtblite.*" $LIB_PATHS
    ;;

  *)
    echo "==================== Linking TBLITE to user paths ===================="
    pkg_install_dir="$with_tblite"
    check_dir "${pkg_install_dir}/include"
    ;;

esac

if [ "$with_tblite" != "__DONTUSE__" ]; then

  TBLITE_DFLAGS="-D__TBLITE -D__DFTD4"
  TBLITE_LIBS="-ltblite -ldftd4 -ls-dftd3 -lmulticharge -lmctc-lib -ltoml-f"

  cat << EOF > "${BUILDDIR}/setup_tblite"
export TBLITE_VER="${tblite_ver}"
EOF

  if [ "$with_tblite" != "__SYSTEM__" ]; then
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "tomlf.mod")
    TOMLF=${TEMP_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "multicharge.mod")
    MCHARGE=${TEMP_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "mstore.mod")
    MSTORE=${TEMP_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "mctc_io.mod")
    MCTC=${TEMP_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "dftd3.mod")
    SDFTD3=${TEMP_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "dftd4.mod")
    DFTD4=${TEMP_LOC%/*}
    TEMP_LOC=$(find ${pkg_install_dir}/include -name "tblite_xtb.mod")
    TBLITE=${TEMP_LOC%/*}

    TBLITE_INCLUDE_DIRS="${pkg_install_dir}/include"
    TBLITE_LINK_LIBRARIES="${pkg_install_dir}/lib"
    TBLITE_CFLAGS="-I'${TOMLF}' -I'${MCTC}' -I'${SDFTD3}' -I'${DFTD4}' -I'${TBLITE}'"
    TBLITE_LDFLAGS="-L'${TBLITE_LINK_LIBRARIES}' -Wl,-rpath,'${TBLITE_LINK_LIBRARIES}'"

    cat << EOF >> "${BUILDDIR}/setup_tblite"
prepend_path LD_LIBRARY_PATH "${TBLITE_LINK_LIBRARIES}"
prepend_path LD_RUN_PATH "${TBLITE_LINK_LIBRARIES}"
prepend_path LIBRARY_PATH "${TBLITE_LINK_LIBRARIES}"
prepend_path CPATH "${TBLITE_INCLUDE_DIRS}"
prepend_path PKG_CONFIG_PATH "${TBLITE_LINK_LIBRARIES}/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi

  cat << EOF >> "${BUILDDIR}/setup_tblite"
export TOMLF="${TOMLF}"
export MCHARGE="${MCHARGE}"
export MSTORE="${MSTORE}"
export MCTC="${MCTC}"
export SDFTD3="${SDFTD3}"
export DFTD4="${DFTD4}"
export TBLITE="${TBLITE}"
export TBLITE_INCLUDE_DIRS="${TBLITE_INCLUDE_DIRS}"
export TBLITE_LINK_LIBRARIES="${TBLITE_LINK_LIBRARIES}"
export TBLITE_ROOT="${pkg_install_dir}"
export TBLITE_DFLAGS="${TBLITE_DFLAGS}" 
export TBLITE_CFLAGS="${TBLITE_CFLAGS}"
export TBLITE_LDFLAGS="${TBLITE_LDFLAGS}"
export TBLITE_LIBS="${TBLITE_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} \${TBLITE_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} \${TBLITE_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} \${TBLITE_LDFLAGS}"
export CP_LIBS="\${TBLITE_LIBS} \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_tblite" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_tblite"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "tblite"
