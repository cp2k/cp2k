#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

cosma_ver="2.8.5"
cosma_sha256="0ad8e42807cd52ca0e92d27304e2a178e449c88cbaf0d94a410e2487369a3fe1"
costa_ver="2.3.3"
costa_sha256="e007ec640c988fb04f5659228f7cc7733bb3c4621cd1fa65f8671420b6807852"
tiled_mm_ver="2.3.2"
tiled_mm_sha256="1f91ca02f6ee8e400835fa90630618baf86a7b425b4bbbb4151068f72658b858"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_cosma" ] && rm "${BUILDDIR}/setup_cosma"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

# HIP-build of Tiled-MM and COSMA requires rocblas
if [ "${with_cosma}" != "__DONTUSE__" ] && [ "${ENABLE_HIP}" = "__TRUE__" ]; then
  check_lib -lrocblas "rocm"
fi

case "${with_cosma}" in
  __INSTALL__)
    require_env OPENBLAS_ROOT
    require_env SCALAPACK_ROOT

    echo "==================== Installing COSMA ===================="
    pkg_install_dir="${INSTALLDIR}/COSMA-${cosma_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "COSMA-${cosma_ver} is already installed, skipping it."
    else
      retrieve_package "${cosma_sha256}" "COSMA-v${cosma_ver}.tar.gz"
      retrieve_package "${costa_sha256}" "COSTA-v${costa_ver}.tar.gz"
      [ -d COSMA-${cosma_ver} ] && rm -rf COSMA-${cosma_ver}
      tar -xzf COSMA-v${cosma_ver}.tar.gz
      [ -d COSTA-${costa_ver} ] && rm -rf COSTA-${costa_ver}
      tar -xzf COSTA-v${costa_ver}.tar.gz
      if [ "${ENABLE_HIP}" = "__TRUE__" ] || [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        retrieve_package "${tiled_mm_sha256}" "Tiled-MM-v${tiled_mm_ver}.tar.gz"
        [ -d Tiled-MM-${tiled_mm_ver} ] && rm -rf Tiled-MM-${tiled_mm_ver}
        tar -xzf Tiled-MM-v${tiled_mm_ver}.tar.gz
      fi
      echo "Installing from scratch into ${pkg_install_dir}"

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
      cmake -B build-cpu \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DBUILD_SHARED_LIBS=NO \
        -DCOSTA_SCALAPACK=${cosma_sl} \
        > cmake-cpu.log 2>&1 || tail_excerpt cmake-cpu.log
      cmake --build build-cpu -t install -j $(get_nprocs) \
        > build-cpu.log 2>&1 || tail_excerpt build-cpu.log

      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        cmake -B build-cuda \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSTA_SCALAPACK=${cosma_sl} \
          > cmake-cuda.log 2>&1 || tail_excerpt cmake-cuda.log
        cmake --build build-cuda -t install -j $(get_nprocs) \
          > build-cuda.log 2>&1 || tail_excerpt build-cuda.log
      fi

      if [ "${ENABLE_HIP}" = "__TRUE__" ]; then
        cmake -B build-hip \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSTA_BLAS=${cosma_blas} \
          -DCOSTA_SCALAPACK=${cosma_sl} \
          > cmake-hip.log 2>&1 || tail_excerpt cmake-hip.log
        cmake --build build-hip -t install -j $(get_nprocs) \
          > build-hip.log 2>&1 || tail_excerpt build-hip.log
      fi
      cd ..

      if [ "${ENABLE_HIP}" = "__TRUE__" ] || [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        cd Tiled-MM-${tiled_mm_ver}
        if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
          cmake -B build-cuda \
            -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
            -DCMAKE_INSTALL_LIBDIR=lib \
            -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
            -DBUILD_SHARED_LIBS=NO \
            -DTILEDMM_GPU_BACKEND=CUDA \
            -DTILEDMM_WITH_EXAMPLES=OFF \
            -DTILEDMM_WITH_TESTS=OFF \
            > cmake-cuda.log 2>&1 || tail_excerpt cmake-cuda.log
          cmake --build build-cuda -t install -j $(get_nprocs) \
            > build-cuda.log 2>&1 || tail_excerpt build-cuda.log
        fi

        if [ "${ENABLE_HIP}" = "__TRUE__" ]; then
          mkdir build-hip && cd build-hip
          cmake -B build-hip \
            -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
            -DCMAKE_INSTALL_LIBDIR=lib \
            -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
            -DBUILD_SHARED_LIBS=NO \
            -DTILEDMM_GPU_BACKEND=ROCM \
            -DTILEDMM_WITH_EXAMPLES=OFF \
            -DTILEDMM_WITH_TESTS=OFF \
            > cmake-hip.log 2>&1 || tail_excerpt cmake-hip.log
          cmake --build build-hip -t install -j $(get_nprocs) \
            > build-hip.log 2>&1 || tail_excerpt build-hip.log
        fi
        cd ..
      fi

      cd COSMA-${cosma_ver}
      cmake -B build-cpu \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DBUILD_SHARED_LIBS=NO \
        -DCOSMA_BLAS=${cosma_blas} \
        -DCOSMA_SCALAPACK=${cosma_sl} \
        -DCOSMA_WITH_TESTS=NO \
        -DCOSMA_WITH_APPS=NO \
        -DCOSMA_WITH_BENCHMARKS=NO \
        > cmake-cpu.log 2>&1 || tail_excerpt cmake-cpu.log
      cmake --build build-cpu -t install -j $(get_nprocs) \
        > build-cpu.log 2>&1 || tail_excerpt build-cpu.log

      # Build CUDA version.
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        cmake -B build-cuda \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSMA_BLAS=CUDA \
          -DCOSMA_SCALAPACK=${cosma_sl} \
          -DCOSMA_WITH_TESTS=NO \
          -DCOSMA_WITH_APPS=NO \
          -DCOSMA_WITH_BENCHMARKS=NO \
          > cmake-cuda.log 2>&1 || tail_excerpt cmake-cuda.log
        cmake --build build-cuda -t install -j $(get_nprocs) \
          > build-cuda.log 2>&1 || tail_excerpt build-cuda.log
      fi

      # Build HIP version.
      if [ "${ENABLE_HIP}" = "__TRUE__" ]; then
        cmake -B build-hip \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DBUILD_SHARED_LIBS=NO \
          -DCOSMA_BLAS=ROCM \
          -DCOSMA_SCALAPACK=${cosma_sl} \
          -DCOSMA_WITH_TESTS=NO \
          -DCOSMA_WITH_APPS=NO \
          -DCOSMA_WITH_BENCHMARKS=NO \
          > cmake-hip.log 2>&1 || tail_excerpt cmake-hip.log
        cmake --build build-hip -t install -j $(get_nprocs) \
          > build-hip.log 2>&1 || tail_excerpt build-hip.log
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding COSMA from system paths ===================="
    check_pkgconfig cosma
    pkg_install_dir="$(pkg-config --variable=prefix cosma)"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking COSMA to user paths ===================="
    pkg_install_dir="${with_cosma}"
    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    COSMA_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && COSMA_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${COSMA_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "${with_cosma}" != "__DONTUSE__" ]; then
  if [ "${with_cosma}" != "__SYSTEM__" ]; then
    pkg_install_dir1="${pkg_install_dir}"
    if [ "${with_cosma}" = "__INSTALL__" ]; then
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        pkg_install_dir1="${pkg_install_dir}-cuda"
      elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
        pkg_install_dir1="${pkg_install_dir}-hip"
      fi
    fi
    cat << EOF > "${BUILDDIR}/setup_cosma"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir1}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir1}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir1}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir1}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir1}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_cosma"
export COSMA_VER="${cosma_ver}"
EOF
  filter_setup "${BUILDDIR}/setup_cosma" "${SETUPFILE}"

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
