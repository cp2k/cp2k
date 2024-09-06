#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

openblas_ver="0.3.28" # Keep in sync with get_openblas_arch.sh
openblas_sha256="f1003466ad074e9b0c8d421a204121100b0751c96fc6fcf3d1456bd12f8a00a1"
openblas_pkg="OpenBLAS-${openblas_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_openblas" ] && rm "${BUILDDIR}/setup_openblas"

OPENBLAS_CFLAGS=""
OPENBLAS_LDFLAGS=""
OPENBLAS_LIBS=""
OPENBLAS_ROOT=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_openblas}" in
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
        download_pkg_from_cp2k_org "${openblas_sha256}" "${openblas_pkg}"
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d OpenBLAS-${openblas_ver} ] && rm -rf OpenBLAS-${openblas_ver}
      tar -zxf ${openblas_pkg}
      cd OpenBLAS-${openblas_ver}

      # Build OpenBLAS with DYNAMIC_ARCH unless native architecture is requested.
      # If the latter fails, then build with DYNAMIC_ARCH. Relying on DYNAMIC_ARCH
      # is more flexible and more reliable on ARM64. However, it increaes the time
      # needed to build OpenBLAS.
      # It is supported to omit the PREFIX when building OpenBLAS as well as
      # omitting build-keys for the install target.
      #
      # NUM_THREADS=128: this is what the most common Linux distros seem to choose atm
      #                  for a good compromise between memory usage and scalability
      # USE_THREADS=1: is not needed since CP2K's ARCH-files do not rely on
      #                plain PThread based OpenBLAS (OpenMP is used).
      #
      # Unfortunately, NO_SHARED=1 breaks ScaLAPACK build.
      BUILD_DYNAMIC=0
      if [ "native" != "${TARGET_CPU}" ] || [ ! "${OPENBLAS_LIBCORE}" ]; then
        BUILD_DYNAMIC=1
      fi
      if [ "0" = "${BUILD_DYNAMIC}" ]; then
        TARGET=$(tr '[:lower:]' '[:upper:]' <<< "${OPENBLAS_LIBCORE}")
        echo "Installing OpenBLAS library for target ${TARGET}"
        if ! make -j $(get_nprocs) \
          TARGET=${TARGET} \
          MAKE_NB_JOBS=0 \
          NUM_THREADS=128 \
          USE_OPENMP=1 \
          NO_AFFINITY=1 \
          CC="${CC}" \
          FC="${FC}" \
          PREFIX="${pkg_install_dir}" \
          > make.${OPENBLAS_LIBCORE}.log 2>&1; then
          tail -n ${LOG_LINES} make.${OPENBLAS_LIBCORE}.log
          BUILD_DYNAMIC=1
        fi
      fi
      if [ "0" != "${BUILD_DYNAMIC}" ]; then
        echo "Installing OpenBLAS library for dynamic target"
        make -j $(get_nprocs) \
          DYNAMIC_ARCH=1 \
          MAKE_NB_JOBS=0 \
          NUM_THREADS=128 \
          USE_OPENMP=1 \
          NO_AFFINITY=1 \
          CC="${CC}" \
          FC="${FC}" \
          PREFIX="${pkg_install_dir}" \
          > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      fi
      make MAKE_NB_JOBS=0 PREFIX="${pkg_install_dir}" install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage2/$(basename ${SCRIPT_NAME})"
    fi
    OPENBLAS_CFLAGS="-I'${pkg_install_dir}/include'"
    OPENBLAS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    OPENBLAS_ROOT="${pkg_install_dir}"
    # Prefer static library if available
    if [ -f "${pkg_install_dir}/lib/libopenblas.a" ]; then
      OPENBLAS_LIBS="-l:libopenblas.a"
    else
      OPENBLAS_LIBS="-lopenblas"
    fi
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
    OPENBLAS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
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
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
prepend_path CPATH "$pkg_install_dir/include"
export OPENBLAS_ROOT=${pkg_install_dir}
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_openblas"
export OPENBLAS_VER="${openblas_ver}"
export OPENBLAS_ROOT="${pkg_install_dir}"
export OPENBLAS_CFLAGS="${OPENBLAS_CFLAGS}"
export OPENBLAS_LDFLAGS="${OPENBLAS_LDFLAGS}"
export OPENBLAS_LIBS="${OPENBLAS_LIBS}"
export MATH_CFLAGS="\${MATH_CFLAGS} ${OPENBLAS_CFLAGS}"
export MATH_LDFLAGS="\${MATH_LDFLAGS} ${OPENBLAS_LDFLAGS}"
export MATH_LIBS="\${MATH_LIBS} ${OPENBLAS_LIBS}"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  cat "${BUILDDIR}/setup_openblas" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_openblas"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "openblas"
