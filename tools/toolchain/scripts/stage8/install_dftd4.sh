#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

dftd4_ver="3.6.0"
dftd4_sha256="0e3e8d5f9e9e5414b9979967c074c953706053832e551d922c27599e7324bace"

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
    echo "==================== Installing GRIMME D4 ===================="
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
      cd dftd4-${dftd4_ver}

      rm -Rf build
      mkdir build
      cd build

      CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${OPENBLAS_ROOT}" cmake \
        -B . -G Ninja \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_C_COMPILER="${MPICC}" \
        -DCMAKE_Fortran_COMPILER="${MPIFC}" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      cmake --build . -j $(get_nprocs) >> build.log 2>&1 || tail -n ${LOG_LINES} build.log
      cmake --install . >> install.log 2>&1 || tail -n ${LOG_LINES} install.log

      cd ..
      echo "==================== Linking Grimme_D4 to user paths ===================="
    fi
    write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    ;;

  __SYSTEM__)
    echo "==================== Finding dftd4 from system paths ===================="
    check_command pkg-config --modversion dftd4
    add_include_from_paths DFTD4_CFLAGS "dftd4.h" $INCLUDE_PATHS
    add_include_from_paths DFTD4_CFLAGS "dftd4.mod" $INCLUDE_PATHS
    add_include_from_paths DFTD4_CFLAGS "mctc_io.mod" $INCLUDE_PATHS
    add_lib_from_paths DFTD4_LDFLAGS "libdftd4.*" $LIB_PATHS
    ;;

  *)
    echo "==================== Linking dftd4 to user paths ===================="
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
    DFTD4_LOC=$(find ${pkg_install_dir}/include -name "multicharge.mod")
    DFTD4_MCHARGE=${DFTD4_LOC%/*}
    DFTD4_LOC=$(find ${pkg_install_dir}/include -name "mstore.mod")
    DFTD4_STORE=${DFTD4_LOC%/*}
    DFTD4_LOC=$(find ${pkg_install_dir}/include -name "mctc_io.mod")
    DFTD4_MCTC=${DFTD4_LOC%/*}
    DFTD4_LOC=$(find ${pkg_install_dir}/include -name "dftd4.mod")
    DFTD4_DFTD4=${DFTD4_LOC%/*}
    # use the lib64 directory if present
    DFTD4_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && DFTD4_LIBDIR="${pkg_install_dir}/lib64"

    DFTD4_CFLAGS="-I'${pkg_install_dir}/include' -I'${DFTD4_DFTD4}' -I'${DFTD4_MCTC}'"
    DFTD4_LDFLAGS="-L'${DFTD4_LIBDIR}' -Wl,-rpath,'${DFTD4_LIBDIR}'"

    cat << EOF >> "${BUILDDIR}/setup_dftd4"
prepend_path LD_LIBRARY_PATH "${DFTD4_LIBDIR}"
prepend_path LD_RUN_PATH "${DFTD4_LIBDIR}"
prepend_path LIBRARY_PATH "${DFTD4_LIBDIR}"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "${DFTD4_LIBDIR}/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi

  cat << EOF >> "${BUILDDIR}/setup_dftd4"
export DFTD4_DFTD4="${DFTD4_DFTD4}"
export DFTD4_MCTC="${DFTD4_MCTC}"
export DFTD4_LIBDIR="${DFTD4_LIBDIR}"
export DFTD4_INCLUDE_DIR="$pkg_install_dir/include"
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
