#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

cosma_ver="2.6.6"
cosma_sha256="1604be101e77192fbcc5551236bc87888d336e402f5409bbdd9dea900401cc37"
costa_ver="2.2.2"
costa_sha256="e87bc37aad14ac0c5922237be5d5390145c9ac6aef0350ed17d86cb2d994e67c"
tiled_mm_ver="2.3"
tiled_mm_sha256="504c6201f5a9be9741c55036bf8e2656ae3f4bc19996295b264ee5e303c9253c"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_cosma" ] && rm "${BUILDDIR}/setup_cosma"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_cosma" in
  __INSTALL__)
    require_env OPENBLAS_ROOT
    require_env SCALAPACK_ROOT

    echo "==================== Installing COSMA ===================="
    pkg_install_dir="${INSTALLDIR}/COSMA-${cosma_ver}"
    install_lock_file="$pkg_install_dir/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "COSMA-${cosma_ver} is already installed, skipping it."
    else
      if [ -f COSMA-v${cosma_ver}.tar.gz ]; then
        echo "COSMA-v${cosma_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${cosma_sha256}" "COSMA-v${cosma_ver}.tar.gz"
        download_pkg_from_cp2k_org "${costa_sha256}" "COSTA-v${costa_ver}.tar.gz"
        download_pkg_from_cp2k_org "${tiled_mm_sha256}" "Tiled-MM-v${tiled_mm_ver}.tar.gz"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d COSMA-${cosma_ver} ] && rm -rf COSMA-${cosma_ver}
      tar -xzf COSMA-v${cosma_ver}.tar.gz
      [ -d Tiled-MM-${tiled_mm_ver} ] && rm -rf Tiled-MM-${tiled_mm_ver}
      tar -xzf Tiled-MM-v${tiled_mm_ver}.tar.gz
      [ -d COSTA-${costa_ver} ] && rm -rf COSTA-${costa_ver}
      tar -xzf COSTA-v${costa_ver}.tar.gz

      case "$MATH_MODE" in
        mkl)
          cosma_blas="MKL"
          cosma_sl="MKL"
          ;;
        cray)
          cosma_blas="CRAY_LIBSCI"
          cosma_sl="CRAY_LIBSCI"
          ;;
        *)
          cosma_blas="OPENBLAS"
          cosma_sl="CUSTOM"
          ;;
      esac

      cd "COSTA-${costa_ver}"
      [ -d build-cpu ] && rm -Rf build-cpu
      mkdir build-cpu && cd build-cpu
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DBUILD_SHARED_LIBS=NO \
        -DCOSTA_SCALAPACK=${cosma_sl} \
        .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..

      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -Rf build-cuda
        mkdir build-cuda && cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSTA_SCALAPACK=${cosma_sl} \
          .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi

      if [ "$ENABLE_HIP" = "__TRUE__" ]; then
        [ -d build-hip ] && rm -Rf build-hip
        mkdir build-hip && cd build-hip
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSTA_BLAS=${cosma_blas} \
          -DCOSTA_SCALAPACK=${cosma_sl} \
          .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi
      cd ..

      cd Tiled-MM-${tiled_mm_ver}
      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        mkdir build-cuda && cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DTILEDMM_GPU_BACKEND=CUDA \
          -DTILEDMM_WITH_EXAMPLES=OFF \
          -DTILEDMM_WITH_TESTS=OFF \
          .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi

      if [ "$ENABLE_HIP" = "__TRUE__" ]; then
        mkdir build-hip && cd build-hip
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DTILEDMM_GPU_BACKEND=ROCM \
          -DTILEDMM_WITH_EXAMPLES=OFF \
          -DTILEDMM_WITH_TESTS=OFF \
          .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi
      cd ..

      cd COSMA-${cosma_ver} && mkdir build-cpu && cd build-cpu
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DBUILD_SHARED_LIBS=NO \
        -DCOSMA_BLAS=${cosma_blas} \
        -DCOSMA_SCALAPACK=${cosma_sl} \
        -DCOSMA_WITH_TESTS=NO \
        -DCOSMA_WITH_APPS=NO \
        -DCOSMA_WITH_BENCHMARKS=NO .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..

      # Build CUDA version.
      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        mkdir build-cuda
        cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSMA_BLAS=CUDA \
          -DCOSMA_SCALAPACK=${cosma_sl} \
          -DCOSMA_WITH_TESTS=NO \
          -DCOSMA_WITH_APPS=NO \
          -DCOSMA_WITH_BENCHMARKS=NO .. \
          > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi

      # Build HIP version.
      if [ "$ENABLE_HIP" = "__TRUE__" ] && $(check_lib -lrocblas "rocm" &> /dev/null); then
        [ -d build-hip ] && rm -rf "build-hip"
        mkdir build-hip
        cd build-hip
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSMA_BLAS=ROCM \
          -DCOSMA_SCALAPACK=${cosma_sl} \
          -DCOSMA_WITH_TESTS=NO \
          -DCOSMA_WITH_APPS=NO \
          -DCOSMA_WITH_BENCHMARKS=NO .. \
          > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    COSMA_ROOT="${pkg_install_dir}"
    COSMA_CFLAGS="-I'${pkg_install_dir}/include'"

    # Check if COSMA is compiled with 64bits and set up COSMA_LIBDIR accordingly.
    COSMA_LIBDIR="${pkg_install_dir}/lib"
    COSMA_LDFLAGS="-L'${COSMA_LIBDIR}' -Wl,-rpath,'${COSMA_LIBDIR}'"
    COSMA_CUDA_LIBDIR="${pkg_install_dir}-cuda/lib"
    COSMA_CUDA_LDFLAGS="-L'${COSMA_CUDA_LIBDIR}' -Wl,-rpath,'${COSMA_CUDA_LIBDIR}'"
    COSMA_HIP_LIBDIR="${pkg_install_dir}-hip/lib"
    COSMA_HIP_LDFLAGS="-L'${COSMA_HIP_LIBDIR}' -Wl,-rpath,'${COSMA_HIP_LIBDIR}'"
    ;;
  __SYSTEM__)
    echo "==================== Finding COSMA from system paths ===================="
    check_command pkg-config --modversion cosma
    add_include_from_paths COSMA_CFLAGS "cosma.h" $INCLUDE_PATHS
    add_lib_from_paths COSMA_LDFLAGS "libcosma.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking COSMA to user paths ===================="
    pkg_install_dir="$with_cosma"

    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    COSMA_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && COSMA_LIBDIR="${pkg_install_dir}/lib64"

    check_dir "$pkg_install_dir/lib"
    check_dir "$pkg_install_dir/include"

    COSMA_CFLAGS="-I'${pkg_install_dir}/include'"
    COSMA_LDFLAGS="-L'${COSMA_LIBDIR}' -Wl,-rpath,'${COSMA_LIBDIR}'"
    ;;
esac
if [ "$with_cosma" != "__DONTUSE__" ]; then
  COSMA_LIBS="-lcosma_prefixed_pxgemm -lcosma -lcosta IF_CUDA(-lTiled-MM|)"
  if [ "$ENABLE_HIP" = "__TRUE__" ] && $(check_lib -lrocblas "rocm" &> /dev/null); then
    COSMA_LIBS+=" IF_HIP(-lTiled-MM|)"
  fi
  if [ "$with_cosma" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_cosma"
prepend_path LD_LIBRARY_PATH "${COSMA_LIBDIR}"
prepend_path LD_RUN_PATH "${COSMA_LIBDIR}"
prepend_path LIBRARY_PATH "${COSMA_LIBDIR}"
prepend_path CPATH "$pkg_install_dir/include"
export COSMA_INCLUDE_DIR="$pkg_install_dir/include"
export COSMA_ROOT="${pkg_install_dir}"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${COSMA_LIBDIR}/pkgconfig"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_cosma"
export COSMA_CFLAGS="${COSMA_CFLAGS}"
export COSMA_LDFLAGS="${COSMA_LDFLAGS}"
export COSMA_CUDA_LDFLAGS="${COSMA_CUDA_LDFLAGS}"
export COSMA_HIP_LDFLAGS="${COSMA_HIP_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__COSMA|)"
export CP_CFLAGS="\${CP_CFLAGS} ${COSMA_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_CUDA(${COSMA_CUDA_LDFLAGS}|IF_HIP(${COSMA_HIP_LDFLAGS}|${COSMA_LDFLAGS}))"
export COSMA_LIBS="${COSMA_LIBS}"
export COSMA_ROOT="$pkg_install_dir"
export COSMA_INCLUDE_DIR="$pkg_install_dir/include"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${COSMA_LIBDIR}/pkgconfig"
export COSMA_VERSION=${cosma_ver}
export CP_LIBS="IF_MPI(${COSMA_LIBS}|) \${CP_LIBS}"
EOF
  cat "${BUILDDIR}/setup_cosma" >> $SETUPFILE

  cat << EOF >> ${INSTALLDIR}/lsan.supp
# leaks related to COSMA (probably, only the last one is actually needed)
leak:cosma::communicator::communicator
leak:cosma::cosma_context<double>::register_state
leak:cosma::pxgemm<double>
leak:cosma::cosma_context<std::complex<double> >::register_state
EOF
fi

load "${BUILDDIR}/setup_cosma"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "cosma"
