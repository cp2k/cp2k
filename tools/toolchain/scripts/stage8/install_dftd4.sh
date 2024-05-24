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
        wget "https://github.com/dftd4/dftd4/archive/refs/tags/v${dftd4_ver}.tar.gz" -O "dftd4-${dftd4_ver}.tar.gz"
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d dftd4-${dftd4_ver} ] && rm -rf dftd4-${dftd4_ver}
      tar -xzf dftd4-${dftd4_ver}.tar.gz
      cd dftd4-${dftd4_ver}

      rm -Rf build
      mkdir build
      cd build

      if [ -n "${MKL_LIBS}" ]; then
        EXTRA_CMAKE_FLAGS=" -DMKLROOT=${MKLROOT} "
      fi
      
      CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${OPENBLAS_ROOT}" cmake \
        -B . -G Ninja \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_CXX_COMPILER="${MPICXX}" \
        -DCMAKE_C_COMPILER="${MPICC}" \
        -DCMAKE_Fortran_COMPILER="${MPIFC}" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        ${EXTRA_CMAKE_FLAGS} .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log

      cmake --build . >> build.log 2>&1 || tail -n ${LOG_LINES} build.log

      cmake --install .  >> install.log 2>&1

      cd ..
      echo "==================== Linking Grimme_D4 to user paths ===================="
    fi
    write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    ;;

  *)
    echo "__INVALID__"
    ;; 

  esac

if [ "$with_dftd4" != "__DONTUSE__" ]; then
  DFTD4_LOC=$(find ${pkg_install_dir}/include -name "dftd4.mod")
  DFTD4_MOD=${DFTD4_LOC%/*}
  DFTD4_LIBS="-ldftd4"
  DFTD4_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
  DFTD4_CFLAGS="-I'${pkg_install_dir}/include' -I'${DFTD4_MOD}'"

  if [ "$with_dftd4" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_dftd4"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi

  cat << EOF >> "${BUILDDIR}/setup_dftd4"
export DFTD4_LIBS="${DFTD4_LIBS}"
export DFTD4_LDFLAGS="${DFTD4_LDFLAGS}"
export DFTD4_CFLAGS="${DFTD4_CFLAGS}"
export DFTD4_MOD="${DFTD4_MOD}"
export DFTD4_LIBRARIES="${pkg_install_dir}/lib"
export DFTD4_INCLUDE_DIR="$pkg_install_dir/include"
export DFTD4_ROOT="${pkg_install_dir}" 
export CP_DFLAGS="\${CP_DFLAGS} -D__DFTD4"
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
