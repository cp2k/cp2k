#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

tblite_ver="0.6.0"
tblite_sha256="372281aedb89234168d00eb691addb303197a9462a9c55d145c835f2cf5e8b42"
tblite_sdftd3_ver="1.4.0"
tblite_dftd4_ver="4.2.0"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_tblite" ] && rm "${BUILDDIR}/setup_tblite"

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
      #retrieve_package "${tblite_sha256}" "tblite-${tblite_ver}.tar.xz"
      retrieve_github_archive "${tblite_sha256}" "tblite-${tblite_ver}.tar.xz" \
        "https://github.com/tblite/tblite/releases/download/v${tblite_ver}" \
        "tblite-${tblite_ver}.tar.xz"
      [ -d tblite-${tblite_ver} ] && rm -rf tblite-${tblite_ver}
      tar -xJf tblite-${tblite_ver}.tar.xz
      cd tblite-${tblite_ver}

      patch -l -d subprojects/s-dftd3 -p1 < "${SCRIPT_DIR}/stage8/simple-dftd3-${tblite_sdftd3_ver}-gradient-fixes.patch" \
        > simple_dftd3_gradient_fixes.patch.log 2>&1 || tail_excerpt simple_dftd3_gradient_fixes.patch.log
      patch -l -d subprojects/dftd4 -p1 < "${SCRIPT_DIR}/stage8/dftd4-${tblite_dftd4_ver}-gradient-fixes.patch" \
        > dftd4_gradient_fixes.patch.log 2>&1 || tail_excerpt dftd4_gradient_fixes.patch.log

      rm -Rf build
      mkdir build
      cd build

      CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${OPENBLAS_ROOT}" cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_TESTING=OFF \
        -DWITH_TESTS=OFF \
        .. \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
      cd ..
    fi
    ;;

  __SYSTEM__)
    echo "==================== Finding tblite from system paths ===================="
    check_command pkg-config --modversion tblite
    TBLITE_INCLUDE_PATH=$(pkg-config --cflags tblite | awk '{print $1}' | cut -dI -f2)
    pkg_install_dir=$(dirname ${TBLITE_INCLUDE_PATH})
    add_include_from_paths TBLITE_CFLAGS "tblite.h" $TBLITE_INCLUDE_PATH
    add_include_from_paths TBLITE_CFLAGS "tblite.mod" $TBLITE_INCLUDE_PATH
    add_include_from_paths TBLITE_CFLAGS "dftd4.mod" $TBLITE_INCLUDE_PATH
    add_include_from_paths TBLITE_CFLAGS "mctc_io.mod" $TBLITE_INCLUDE_PATH
    add_include_from_paths TBLITE_CFLAGS "mstore.mod" $TBLITE_INCLUDE_PATH
    add_include_from_paths TBLITE_CFLAGS "multicharge.mod" $TBLITE_INCLUDE_PATH
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
export TBLITE_DFTD4_VER="${tblite_dftd4_ver}"
export TBLITE_MULTICHARGE_VER="${tblite_multicharge_ver}"
export TBLITE_SDFTD3_VER="${tblite_sdftd3_ver}"
export TBLITE_MCTC_VER="${tblite_mctc_ver}"
export TBLITE_TOMLF_VER="${tblite_tomlf_ver}"
EOF

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

  if [ "$with_tblite" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_tblite"
prepend_path PATH "${pkg_install_dir}/bin"
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
  filter_setup "${BUILDDIR}/setup_tblite" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_tblite"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "tblite"
