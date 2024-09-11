#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

spla_ver="1.6.1"
spla_sha256="62b51e6ce05c41cfc1c6f6600410f9549a209c50f0331e1db41047f94493e02f"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_SpLA" ] && rm "${BUILDDIR}/setup_SpLA"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_spla}" in
  __INSTALL__)
    echo "==================== Installing spla ===================="
    pkg_install_dir="${INSTALLDIR}/SpLA-${spla_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "SpLA-${spla_ver} is already installed, skipping it."
    else
      if [ -f SpLA-${spla_ver}.tar.gz ]; then
        echo "SpLA-${spla_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${spla_sha256}" "SpLA-${spla_ver}.tar.gz"

      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d SpLA-${spla_ver} ] && rm -rf SpLA-${spla_ver}
      tar -xzf SpLA-${spla_ver}.tar.gz
      cd spla-${spla_ver}
      mkdir -p build-cpu
      cd build-cpu
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DSPLA_FORTRAN=ON \
        -DSPLA_INSTALL=ON \
        -DSPLA_STATIC=ON \
        .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..

      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        mkdir build-cuda
        cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DSPLA_FORTRAN=ON \
          -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
          -DSPLA_INSTALL=ON \
          -DSPLA_STATIC=ON \
          -DSPLA_GPU_BACKEND=CUDA \
          .. \
          > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        install -d ${pkg_install_dir}/lib/cuda
        [ -f src/libspla.a ] && install -m 644 src/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
        [ -f src/libspla.so ] && install -m 644 src/*.so ${pkg_install_dir}/lib/cuda >> install.log 2>&1
      fi

      if [ "$ENABLE_HIP" = "__TRUE__" ]; then

        case "${GPUVER}" in
          K20X | K40 | K80 | P100 | V100 | A100 | A40 | H100)
            [ -d build-cuda ] && rm -rf "build-cuda"
            mkdir build-cuda
            cd build-cuda
            cmake \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
              -DCMAKE_INSTALL_LIBDIR=lib \
              -DCMAKE_VERBOSE_MAKEFILE=ON \
              -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
              -DSPLA_FORTRAN=ON \
              -DSPLA_INSTALL=ON \
              -DSPLA_STATIC=ON \
              -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
              -DSPLA_GPU_BACKEND=CUDA \
              .. \
              > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
            make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
            install -d ${pkg_install_dir}/lib/hip
            [ -f src/libspla.a ] && install -m 644 src/*.a ${pkg_install_dir}/lib/hip >> install.log 2>&1
            [ -f src/libspla.so ] && install -m 644 src/*.so ${pkg_install_dir}/lib/hip >> install.log 2>&1
            ;;
          Mi50 | Mi100 | Mi200 | Mi250)
            [ -d build-hip ] && rm -rf "build-hip"
            mkdir build-hip
            cd build-hip
            cmake \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
              -DCMAKE_INSTALL_LIBDIR=lib \
              -DCMAKE_VERBOSE_MAKEFILE=ON \
              -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
              -DSPLA_FORTRAN=ON \
              -DSPLA_INSTALL=ON \
              -DSPLA_STATIC=ON \
              -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
              -DSPLA_GPU_BACKEND=ROCM \
              .. \
              > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
            make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
            install -d ${pkg_install_dir}/lib/hip
            [ -f src/libspla.a ] && install -m 644 src/*.a ${pkg_install_dir}/lib/hip >> install.log 2>&1
            [ -f src/libspla.so ] && install -m 644 src/*.so ${pkg_install_dir}/lib/hip >> install.log 2>&1
            ;;
          *) ;;
        esac
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    SPLA_ROOT="${pkg_install_dir}"
    SPLA_CFLAGS="-I'${pkg_install_dir}/include/spla'"
    SPLA_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    SPLA_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath,'${pkg_install_dir}/lib/cuda'"
    SPLA_HIP_LDFLAGS="-L'${pkg_install_dir}/lib/hip' -Wl,-rpath,'${pkg_install_dir}/lib/hip'"
    ;;
  __SYSTEM__)
    echo "==================== Finding spla from system paths ===================="
    check_command pkg-config --modversion spla
    add_include_from_paths SPLA_CFLAGS "spla.h" $INCLUDE_PATHS
    add_lib_from_paths SPLA_LDFLAGS "libspla.*" $LIB_PATHS
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking spla to user paths ===================="
    pkg_install_dir="$with_spla"

    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    SPLA_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && SPLA_LIBDIR="${pkg_install_dir}/lib64"

    check_dir "${SPLA_LIBDIR}"
    check_dir "${pkg_install_dir}/include/spla"
    SPLA_CFLAGS="-I'${pkg_install_dir}/include/spla'"
    SPLA_LDFLAGS="-L'${SPLA_LIBDIR}' -Wl,-rpath,'${SPLA_LIBDIR}'"
    ;;
esac
if [ "$with_spla" != "__DONTUSE__" ]; then
  SPLA_LIBS="-lspla"
  if [ "$with_spla" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_spla"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include/spla"
export SPLA_INCLUDE_DIR="$pkg_install_dir/include/spla"
export SPLA_LIBS="-lspla"
export SPLA_ROOT="${pkg_install_dir}"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_spla"
export SPLA_VER="${spla_ver}"
export SPLA_CFLAGS="${SPLA_CFLAGS}"
export SPLA_LDFLAGS="${SPLA_LDFLAGS}"
export SPLA_CUDA_LDFLAGS="${SPLA_CUDA_LDFLAGS}"
export SPLA_HIP_LDFLAGS="${SPLA_HIP_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_HIP(-D__OFFLOAD_GEMM|) IF_CUDA(-D__OFFLOAD_GEMM|) ${OFFLOAD_DFLAGS} IF_MPI(-D__SPLA|)"
export CP_CFLAGS="\${CP_CFLAGS} ${SPLA_CFLAGS}"
export SPLA_LIBRARY="-lspla"
export SPLA_ROOT="${pkg_install_dir}"
export SPLA_INCLUDE_DIR="${pkg_install_dir}/include/spla"
export CP_LIBS="IF_MPI(${SPLA_LIBS}|) \${CP_LIBS}"
EOF
  if [ "$ENABLE_HIP" = "__TRUE__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_spla"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_HIP(${SPLA_HIP_LDFLAGS}|${SPLA_LDFLAGS})"
EOF
  elif [ "$ENABLE_CUDA" = "__TRUE__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_spla"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_CUDA(${SPLA_CUDA_LDFLAGS}|${SPLA_LDFLAGS})"
EOF
  else
    cat << EOF >> "${BUILDDIR}/setup_spla"
export CP_LDFLAGS="\${CP_LDFLAGS} ${SPLA_LDFLAGS}"
EOF
  fi
  cat "${BUILDDIR}/setup_spla" >> $SETUPFILE
fi

load "${BUILDDIR}/setup_spla"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "spla"
