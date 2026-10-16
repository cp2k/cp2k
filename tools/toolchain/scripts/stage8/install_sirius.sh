#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

sirius_ver="7.12.0"
sirius_sha256="25d14a31faf9cd07b19a327aba5f8f024be63cf6d9d033eada4410727bf95ce9"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_sirius" ] && rm "${BUILDDIR}/setup_sirius"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_sirius" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing SIRIUS ===================="
    pkg_install_dir="${INSTALLDIR}/sirius-${sirius_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "sirius-${sirius_ver} is already installed, skipping it."
    else
      retrieve_package "${sirius_sha256}" "SIRIUS-${sirius_ver}.tar.gz"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d sirius-${sirius_ver} ] && rm -rf sirius-${sirius_ver}
      tar -xzf SIRIUS-${sirius_ver}.tar.gz
      cd SIRIUS-${sirius_ver}
      if [ "${with_elpa}" != "__DONTUSE__" ]; then
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_ELPA=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      if [ "${math_mode}" = "mkl" ]; then
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_MKL=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      if [ "${with_tblite}" != "__DONTUSE__" ]; then
        # tblite includes s-dftd3
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_DFTD3=ON -DSIRIUS_USE_DFTD4=ON ${EXTRA_CMAKE_FLAGS}"
      elif [ "${with_dftd4}" != "__DONTUSE__" ]; then
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_DFTD4=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      cmake -B build \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_BUILD_TYPE="Release" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DSIRIUS_USE_SCALAPACK=ON \
        -DSIRIUS_USE_VCSQNM=ON \
        -DSIRIUS_USE_VDWXC=ON \
        -DSIRIUS_USE_PUGIXML=ON \
        -DSIRIUS_USE_MEMORY_POOL=OFF \
        ${EXTRA_CMAKE_FLAGS} \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      cmake --build build -t install -j $(get_nprocs) \
        > build.log 2>&1 || tail_excerpt build.log

      # now do we have cuda as well
      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        echo "Installing from scratch into ${pkg_install_dir}/cuda"
        cmake -B build-cuda \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir}/cuda \
          -DCMAKE_INSTALL_LIBDIR="lib" \
          -DCMAKE_BUILD_TYPE="Release" \
          -DCMAKE_CUDA_FLAGS="-allow-unsupported-compiler" \
          -DSIRIUS_USE_CUDA=ON \
          -DCMAKE_CUDA_ARCHITECTURES="${ARCH_NUM}" \
          -DSIRIUS_USE_MEMORY_POOL=OFF \
          -DBUILD_SHARED_LIBS=OFF \
          -DSIRIUS_USE_SCALAPACK=ON \
          -DSIRIUS_USE_VCSQNM=ON \
          -DSIRIUS_USE_PUGIXML=ON \
          -DSIRIUS_USE_VDWXC=ON \
          ${EXTRA_CMAKE_FLAGS} \
          > cmake-cuda.log 2>&1 || tail_excerpt cmake-cuda.log
        cmake --build build-cuda --target install -j $(get_nprocs) \
          > build-cuda.log 2>&1 || tail_excerpt build-cuda.log
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    check_lib -lsirius "sirius"
    check_lib -lsirius_cxx "sirius_cxx"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libsirius.*" $LIB_PATHS)))"
    ;;
  *)
    echo "==================== Linking SIRIUS to user paths ===================="
    pkg_install_dir="${with_sirius}"
    SIRIUS_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && SIRIUS_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${SIRIUS_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_sirius" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_sirius"
export SIRIUS_VER="${sirius_ver}"
EOF
  if [ "$with_sirius" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_sirius"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
    if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
      cat << EOF >> "${BUILDDIR}/setup_sirius"
prepend_path PATH "${pkg_install_dir}/cuda/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/cuda/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/cuda/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/cuda/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}/cuda"
EOF
    fi
  fi
  filter_setup "${BUILDDIR}/setup_sirius" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_sirius"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "sirius"
