#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

openblas_ver="0.3.19" # Keep in sync with get_openblas_arch.sh
openblas_sha256="947f51bfe50c2a0749304fbe373e00e7637600b0a47b78a51382aeb30ca08562"
openblas_pkg="OpenBLAS-${openblas_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_openblas" ] && rm "${BUILDDIR}/setup_openblas"

OPENBLAS_CFLAGS=''
OPENBLAS_LDFLAGS=''
OPENBLAS_LIBS=''
OPENBLASROOT=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_openblas" in
  __INSTALL__)
    echo "==================== Installing OpenBLAS ===================="
    pkg_install_dir="${INSTALLDIR}/openblas-${openblas_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "openblas-${openblas_ver} is already installed, skipping it."
    else
      if [ -f ${openblas_pkg} ]; then
        echo "${openblas_pkg} is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${openblas_sha256} \
          https://www.cp2k.org/static/downloads/${openblas_pkg}
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d OpenBLAS-${openblas_ver} ] && rm -rf OpenBLAS-${openblas_ver}
      tar -zxf ${openblas_pkg}
      cd OpenBLAS-${openblas_ver}

      # First attempt to make openblas using auto detected
      # TARGET, if this fails, then make with forced
      # TARGET=NEHALEM
      #
      # wrt NUM_THREADS=64: this is what the most common Linux distros seem to choose atm
      #                     for a good compromise between memory usage and scalability
      #
      # Unfortunately, NO_SHARED=1 breaks ScaLAPACK build.
      if [ "${generic}" == "__TRUE__" ]; then
        make -j $(get_nprocs) \
          MAKE_NB_JOBS=0 \
          TARGET=NEHALEM \
          NUM_THREADS=64 \
          USE_THREAD=1 \
          USE_OPENMP=1 \
          CC="${CC}" \
          FC="${FC}" \
          PREFIX="${pkg_install_dir}" \
          > make.generic.log 2>&1
      else
        (
          make -j $(get_nprocs) \
            MAKE_NB_JOBS=0 \
            NUM_THREADS=64 \
            USE_THREAD=1 \
            USE_OPENMP=1 \
            CC="${CC}" \
            FC="${FC}" \
            PREFIX="${pkg_install_dir}" \
            > make.log 2>&1
        ) || (
          make -j $(get_nprocs) \
            MAKE_NB_JOBS=0 \
            TARGET=NEHALEM \
            NUM_THREADS=64 \
            USE_THREAD=1 \
            USE_OPENMP=1 \
            CC="${CC}" \
            FC="${FC}" \
            PREFIX="${pkg_install_dir}" \
            > make.nehalem.log 2>&1
        )
      fi
      make -j $(get_nprocs) \
        MAKE_NB_JOBS=0 \
        NUM_THREADS=64 \
        USE_THREAD=1 \
        USE_OPENMP=1 \
        CC="${CC}" \
        FC="${FC}" \
        PREFIX="${pkg_install_dir}" \
        install > install.log 2>&1
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage2/$(basename ${SCRIPT_NAME})"
    fi
    OPENBLAS_CFLAGS="-I'${pkg_install_dir}/include'"
    OPENBLAS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    OPENBLASROOT='${pkg_install_dir}'
    OPENBLAS_LIBS="-lopenblas"
    ;;
  __SYSTEM__)
    echo "==================== Finding LAPACK from system paths ===================="
    # assume that system openblas is threaded
    check_lib -lopenblas "OpenBLAS"
    OPENBLAS_LIBS="-lopenblas"
    # detect separate omp builds
    check_lib -lopenblas_openmp 2> /dev/null && OPENBLAS_LIBS="-lopenblas_openmp"
    check_lib -lopenblas_omp 2> /dev/null && OPENBLAS_LIBS="-lopenblas_omp"
    add_include_from_paths OPENBLAS_CFLAGS "openblas_config.h" $INCLUDE_PATHS
    add_lib_from_paths OPENBLAS_LDFLAGS "libopenblas.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking LAPACK to user paths ===================="
    pkg_install_dir="$with_openblas"
    check_dir "${pkg_install_dir}/include"
    check_dir "${pkg_install_dir}/lib"
    OPENBLAS_CFLAGS="-I'${pkg_install_dir}/include'"
    OPENBLAS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    OPENBLAS_LIBS="-lopenblas"
    # detect separate omp builds
    (__libdir="${pkg_install_dir}/lib" LIB_PATHS="__libdir" check_lib -lopenblas_openmp 2> /dev/null) &&
      OPENBLAS_LIBS="-lopenblas_openmp"
    (__libdir="${pkg_install_dir}/lib" LIB_PATHS="__libdir" check_lib -lopenblas_omp 2> /dev/null) &&
      OPENBLAS_LIBS="-lopenblas_omp"
    ;;
esac
if [ "$with_openblas" != "__DONTUSE__" ]; then
  if [ "$with_openblas" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_openblas"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
export OPENBLASROOT=${pkg_install_dir}
EOF
    cat "${BUILDDIR}/setup_openblas" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_openblas"
export OPENBLASROOT="${pkg_install_dir}"
export OPENBLAS_CFLAGS="${OPENBLAS_CFLAGS}"
export OPENBLAS_LDFLAGS="${OPENBLAS_LDFLAGS}"
export OPENBLAS_LIBS="${OPENBLAS_LIBS}"
export MATH_CFLAGS="\${MATH_CFLAGS} ${OPENBLAS_CFLAGS}"
export MATH_LDFLAGS="\${MATH_LDFLAGS} ${OPENBLAS_LDFLAGS}"
export MATH_LIBS="\${MATH_LIBS} ${OPENBLAS_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_openblas"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "openblas"
