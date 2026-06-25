#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

sirius_ver="7.11.1"
sirius_sha256="ce77dd168b2c3ef4a89cc2d6b163cb00b3d8b9d4b0652bc9de187dc7e0e74d77"

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
    ARCH=$(uname -m)
    SIRIUS_OPT="-O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
    if [ "$ARCH" = "ppc64le" ]; then
      SIRIUS_OPT="-O3 -DNDEBUG -mcpu=power8 -mtune=power8 -funroll-loops -ftree-vectorize  -mvsx  -maltivec  -mpopcntd  -mveclibabi=mass -fvect-cost-model -fpeel-loops -mcmodel=medium ${MATH_CFLAGS}"
      SIRIUS_DBG="-O2 -g -mcpu=power8 -mtune=power8 -funroll-loops -ftree-vectorize  -mvsx  -maltivec  -mpopcntd  -mveclibabi=mass -fvect-cost-model -fpeel-loops -mcmodel=medium ${MATH_CFLAGS}"
    fi

    if [ "$ARCH" = "x86_64" ]; then
      if [ "${with_intel}" != "__DONTUSE__" ]; then
        SIRIUS_OPT="-DNDEBUG -O2 -g ${MATH_CFLAGS}"
        SIRIUS_DBG="-O1 -g ${MATH_CFLAGS}"
        # SIRIUS_DBG and SIRIUS_OPT are not really considered by CMake and rather the CMAKE_BUILD_TYPE matters.
        # The CMAKE_BUILD_TYPEs "Release" and "RelWithDebInfo" employ -O3/-O2, but already -O2 makes the SIRIUS
        # build quite memory and time intensive. The CMAKE_BUILD_TYPE "Debug" allows for fast compilation, but it
        # generates very slow code.
        # EXTRA_CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS= ${EXTRA_CMAKE_FLAGS}"
        EXTRA_CMAKE_FLAGS="-DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_FLAGS= ${EXTRA_CMAKE_FLAGS}"
      else
        SIRIUS_OPT="-O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
        SIRIUS_DBG="-O2 -g -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
      fi
    fi

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

      patch -l -p1 < "${SCRIPT_DIR}/stage8/sirius-linking.patch" \
        > sirius-linking.patch.log 2>&1 || tail_excerpt sirius-linking.patch.log

      rm -Rf build
      mkdir build
      cd build
      if [ "${with_elpa}" != "__DONTUSE__" ]; then
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_ELPA=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      if [ "${math_mode}" == "mkl" ]; then
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_MKL=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      if [ "${with_tblite}" != "__DONTUSE__" ]; then
        # tblite includes s-dftd3
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_DFTD3=ON -DSIRIUS_USE_DFTD4=ON ${EXTRA_CMAKE_FLAGS}"
      elif [ "${with_dftd4}" != "__DONTUSE__" ]; then
        EXTRA_CMAKE_FLAGS="-DSIRIUS_USE_DFTD4=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_CXX_FLAGS_RELEASE="${SIRIUS_OPT}" \
        -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
        -DCMAKE_CXX_COMPILER="${MPICXX}" \
        -DCMAKE_C_COMPILER="${MPICC}" \
        -DCMAKE_Fortran_COMPILER="${MPIFC}" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DSIRIUS_USE_SCALAPACK=ON \
        -DSIRIUS_USE_VCSQNM=ON \
        -DSIRIUS_USE_VDWXC=ON \
        -DSIRIUS_USE_PUGIXML=ON \
        -DSIRIUS_USE_MEMORY_POOL=OFF \
        ${EXTRA_CMAKE_FLAGS} .. \
        > cmake.log 2>&1 || tail_excerpt cmake.log

      make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
      cd ..

      # now do we have cuda as well

      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        mkdir build-cuda
        cd build-cuda
        cmake \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir}/cuda \
          -DCMAKE_INSTALL_LIBDIR="lib" \
          -DCMAKE_CXX_FLAGS_RELEASE="${SIRIUS_OPT}" \
          -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
          -DCMAKE_CUDA_FLAGS="-std=c++14 -allow-unsupported-compiler" \
          -DSIRIUS_USE_CUDA=ON \
          -DGPU_MODEL=P100 \
          -DSIRIUS_USE_MEMORY_POOL=OFF \
          -DBUILD_SHARED_LIBS=OFF \
          -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
          -DSIRIUS_USE_SCALAPACK=ON \
          -DSIRIUS_USE_VCSQNM=ON \
          -DSIRIUS_USE_PUGIXML=ON \
          -DSIRIUS_USE_VDWXC=ON \
          -DCMAKE_CXX_COMPILER="${MPICXX}" \
          -DCMAKE_C_COMPILER="${MPICC}" \
          -DCMAKE_Fortran_COMPILER="${MPIFC}" \
          ${EXTRA_CMAKE_FLAGS} .. \
          >> cmake.log 2>&1 || tail_excerpt cmake.log
        make -j $(get_nprocs) install >> make.log 2>&1 || tail_excerpt make.log
        cd ..
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})" \
        "${SCRIPT_DIR}/stage8/sirius-linking.patch"
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
export SIRIUS_ROOT="${pkg_install_dir}"
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
EOF
    fi
  fi
  filter_setup "${BUILDDIR}/setup_sirius" "${SETUPFILE}"
  cat << EOF >> ${INSTALLDIR}/lsan.supp
# leaks related to SIRIUS
leak:cublasXtDeviceSelect
leak:sirius::sirius_free_object_handler
leak:sirius::sddk::memory_pool::free
leak:sirius::sddk::memory_block_descriptor::free_subblock
EOF
fi

load "${BUILDDIR}/setup_sirius"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "sirius"
