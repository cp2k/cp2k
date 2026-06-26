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
    echo "==================== Installing SpLA ===================="
    pkg_install_dir="${INSTALLDIR}/SpLA-${spla_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "SpLA-${spla_ver} is already installed, skipping it."
    else
      retrieve_package "${spla_sha256}" "SpLA-${spla_ver}.tar.gz"
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
        > cmake.log 2>&1 || tail_excerpt cmake.log
      make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
      cd ..

      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        mkdir build-cuda
        cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-cuda" \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_VERBOSE_MAKEFILE=ON \
          -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
          -DSPLA_FORTRAN=ON \
          -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
          -DSPLA_INSTALL=ON \
          -DSPLA_STATIC=ON \
          -DSPLA_GPU_BACKEND=CUDA \
          .. \
          > cmake.log 2>&1 || tail_excerpt cmake.log
        make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
      fi

      if [ "$ENABLE_HIP" = "__TRUE__" ]; then

        case "${GPUVER}" in
          K20X | K40 | K80 | P100 | V100 | A100 | A40 | H100)
            [ -d build-cuda ] && rm -rf "build-cuda"
            mkdir build-cuda
            cd build-cuda
            cmake \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
              -DCMAKE_INSTALL_LIBDIR=lib \
              -DCMAKE_VERBOSE_MAKEFILE=ON \
              -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
              -DSPLA_FORTRAN=ON \
              -DSPLA_INSTALL=ON \
              -DSPLA_STATIC=ON \
              -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
              -DSPLA_GPU_BACKEND=CUDA \
              .. \
              > cmake.log 2>&1 || tail_excerpt cmake.log
            make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
            ;;
          Mi50 | Mi100 | Mi200 | Mi250)
            [ -d build-hip ] && rm -rf "build-hip"
            mkdir build-hip
            cd build-hip
            cmake \
              -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}-hip" \
              -DCMAKE_INSTALL_LIBDIR=lib \
              -DCMAKE_VERBOSE_MAKEFILE=ON \
              -DCMAKE_BUILD_TYPE="RelWithDebInfo" \
              -DSPLA_FORTRAN=ON \
              -DSPLA_INSTALL=ON \
              -DSPLA_STATIC=ON \
              -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
              -DSPLA_GPU_BACKEND=ROCM \
              .. \
              > cmake.log 2>&1 || tail_excerpt cmake.log
            make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
            ;;
          *) ;;
        esac
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding SpLA from system paths ===================="
    check_pkgconfig spla
    pkg_install_dir="$(pkg-config --variable=prefix spla)"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking SpLA to user paths ===================="
    pkg_install_dir="$with_spla"
    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    SPLA_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && SPLA_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${SPLA_LIBDIR}"
    check_dir "${pkg_install_dir}/include/spla"
    ;;
esac
if [ "$with_spla" != "__DONTUSE__" ]; then
  if [ "$with_spla" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_spla"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_spla"
export SPLA_VER="${spla_ver}"
export SPLA_ROOT="${pkg_install_dir}"
EOF
  filter_setup "${BUILDDIR}/setup_spla" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_spla"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "spla"
