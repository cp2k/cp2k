#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

sirius_ver="7.3.1"
sirius_sha256="8bf9848b8ebf0b43797fd359adf8c84f00822de4eb677e3049f22baa72735e98"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

if [ "$MPI_MODE" = "no" ]; then
  report_warning $LINENO "MPI is disabled, skipping sirius installation"
  echo 'with_sirius="__FALSE__"' >> ${BUILDDIR}/setup_sirius
  exit 0
fi

[ -f "${BUILDDIR}/setup_sirius" ] && rm "${BUILDDIR}/setup_sirius"

SIRIUS_CFLAGS=''
SIRIUS_LDFLAGS=''
SIRIUS_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_sirius" in
  __DONTUSE__) ;;

  __INSTALL__)
    echo "==================== Installing SIRIUS ===================="
    require_env FFTW_LDFLAGS
    require_env FFTW_LIBS
    require_env FFTW_CFLAGS
    require_env ELPA_ROOT
    require_env ELPA_LDFLAGS
    require_env ELPA_LIBS
    require_env ELPA_CFLAGS
    require_env GSL_ROOT
    require_env GSL_LDFLAGS
    require_env GSL_CFLAGS
    require_env GSL_LIBS
    require_env GSL_INCLUDE_DIR
    require_env GSL_LIBRARY
    require_env GSL_CBLAS_LIBRARY
    require_env MATH_LIBS
    require_env MPI_LDFLAGS
    require_env MPI_LIBS
    require_env SCALAPACK_ROOT
    require_env SCALAPACK_LDFLAGS
    require_env SCALAPACK_CFLAGS
    require_env SCALAPACK_LIBS
    require_env LIBXC_LIBS
    require_env LIBXC_CFLAGS
    require_env LIBXC_LDFLAGS
    require_env SPGLIB_LIBS
    require_env SPGLIB_CFLAGS
    require_env SPGLIB_LDFLAGS
    require_env HDF5_LIBS
    require_env HDF5_CFLAGS
    require_env HDF5_LDFLAGS
    require_env LIBVDWXC_CFLAGS
    require_env LIBVDWXC_LIBS
    require_env LIBVDWXC_LDFLAGS
    require_env SPFFT_ROOT
    require_env SPFFT_CFLAGS
    require_env SPFFT_LDFLAGS
    require_env SPFFT_LIBS
    require_env SPLA_ROOT
    require_env SPLA_CFLAGS
    require_env SPLA_LDFLAGS
    require_env SPLA_LIBS
    ARCH=$(uname -m)
    SIRIUS_OPT="-O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
    if [ "$ARCH" = "ppc64le" ]; then
      SIRIUS_OPT="-O3 -DNDEBUG -mcpu=power8 -mtune=power8 -funroll-loops -ftree-vectorize  -mvsx  -maltivec  -mpopcntd  -mveclibabi=mass -fvect-cost-model -fpeel-loops -mcmodel=medium ${MATH_CFLAGS}"
      SIRIUS_DBG="-O2 -g -mcpu=power8 -mtune=power8 -funroll-loops -ftree-vectorize  -mvsx  -maltivec  -mpopcntd  -mveclibabi=mass -fvect-cost-model -fpeel-loops -mcmodel=medium ${MATH_CFLAGS}"
    fi

    if [ "$ARCH" = "x86_64" ]; then
      SIRIUS_OPT="-O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
      SIRIUS_DBG="-O2 -g -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}"
    fi

    pkg_install_dir="${INSTALLDIR}/sirius-${sirius_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "sirius_dist-${sirius_ver} is already installed, skipping it."
    else
      if [ -f SIRIUS-${sirius_ver}.tar.gz ]; then
        echo "sirius_${sirius_ver}.tar.gz is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${sirius_sha256} \
          "https://github.com/electronic-structure/SIRIUS/archive/v${sirius_ver}.tar.gz" \
          -o SIRIUS-${sirius_ver}.tar.gz
      fi

      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d sirius-${sirius_ver} ] && rm -rf sirius-${sirius_ver}
      tar -xzf SIRIUS-${sirius_ver}.tar.gz
      cd SIRIUS-${sirius_ver}

      rm -Rf build
      mkdir build
      cd build
      # if [ -n "$ELPA_LIBS" ] ; then
      #     if [ -s "$ELPA_ROOT" ] ; then
      #         export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$ELPA_ROOT/lib/pkgconfig:$ELPA_ROOT/lib64/pkgconfig
      #     fi
      #     EXTRA_CMAKE_FLAGS="-DUSE_ELPA=ON -DELPA_INCLUDE_DIR=${ELPA_ROOT}/include/elpa-${ELPA_VERSION} ${EXTRA_CMAKE_FLAGS}"
      # fi

      if [ -n "${SCALAPACK_LIBS}" ]; then
        export SCALAPACK_LIB="${SCALAPACK_LIBS}"
        if [ -s "${SCALAPACK_ROOT}" ]; then
          EXTRA_CMAKE_FLAGS="-DUSE_SCALAPACK=ON -DSCALAPACK_INCLUDE_DIR=${SCALAPACK_ROOT}/include ${EXTRA_CMAKE_FLAGS}"
        else
          EXTRA_CMAKE_FLAGS="-DUSE_SCALAPACK=ON ${EXTRA_CMAKE_FLAGS}"
        fi
      fi
      if [ -n "${HDF5_LIBS}" ]; then
        CMAKE_PREFIX_PATH="${HDF5_ROOT} ${CMAKE_PREFIX_PATH}"
      fi
      if [ -n "${LIBVDWXC_LIBS}" ]; then
        CMAKE_PREFIX_PATH="${LIBVDWXC_ROOT} ${CMAKE_PREFIX_PATH}"
        EXTRA_CMAKE_FLAGS="-DUSE_VDWXC=ON ${EXTRA_CMAKE_FLAGS}"
      else
        EXTRA_CMAKE_FLAGS="-DUSE_VDWXC=OFF ${EXTRA_CMAKE_FLAGS}"
      fi
      if [ -n "${MKL_LIBS}" ]; then
        EXTRA_CMAKE_FLAGS="-DUSE_MKL=ON -DMKL_DEF_LIBRARY=${MKLROOT}/lib/intel64 -DUSE_SCALAPACK=ON ${EXTRA_CMAKE_FLAGS}"
      fi
      SpFFT_DIR="${SpFFT_ROOT}/lib/cmake/SpFFT"
      SpLA_DIR="${SpLA_ROOT}/lib/cmake/SPLA"
      CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${GSL_ROOT}:${SPGLIB_ROOT}:${LIBXC_ROOT}:${SpFFT_DIR}:${SpLA_DIR}" cmake \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_CXX_FLAGS_RELEASE="${SIRIUS_OPT}" \
        -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
        -DCMAKE_CXX_COMPILER="${MPICXX}" \
        -DCMAKE_C_COMPILER="${MPICC}" \
        -DCMAKE_Fortran_COMPILER="${MPIFC}" \
        -DCMAKE_VERBOSE_MAKEFILE=ON \
        -DBUILD_SHARED_LIBS=OFF \
        -DUSE_ELPA=OFF \
        ${EXTRA_CMAKE_FLAGS} .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log

      make -j $(get_nprocs) -C src >> make.log 2>&1 || tail -n ${LOG_LINES} make.log

      install -d "${pkg_install_dir}/include" >> install.log 2>&1
      install -d "${pkg_install_dir}/lib" >> install.log 2>&1
      cp -R ../src/* "${pkg_install_dir}/include" >> install.log 2>&1
      install -m 644 src/*.a "${pkg_install_dir}/lib" >> install.log 2>&1
      install -m 644 src/mod_files/*.mod "${pkg_install_dir}/include" >> install.log 2>&1
      cd ..

      # now do we have cuda as well

      if [ "$ENABLE_CUDA" = "__TRUE__" ]; then
        [ -d build-cuda ] && rm -rf "build-cuda"
        mkdir build-cuda
        cd build-cuda
        CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}:${GSL_ROOT}:${SPGLIB_ROOT}:${LIBXC_ROOT}:${SpFFT_DIR}:${SpLA_DIR}" cmake \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
          -DCMAKE_CXX_FLAGS_RELEASE="${SIRIUS_OPT}" \
          -DCMAKE_CXX_FLAGS_RELWITHDEBINFO="${SIRIUS_DBG}" \
          -DCMAKE_CUDA_FLAGS="-std=c++14 -allow-unsupported-compiler" \
          -DUSE_CUDA=ON \
          -DUSE_ELPA=OFF \
          -DGPU_MODEL=P100 \
          -DBUILD_SHARED_LIBS=OFF \
          -DCMAKE_CXX_COMPILER="${MPICXX}" \
          -DCMAKE_C_COMPILER="${MPICC}" \
          -DCMAKE_Fortran_COMPILER="${MPIFC}" \
          ${EXTRA_CMAKE_FLAGS} .. \
          >> cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) -C src >> make.log 2>&1 || tail -n ${LOG_LINES} make.log
        install -d ${pkg_install_dir}/lib/cuda
        install -d ${pkg_install_dir}/include/cuda
        install -m 644 src/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
        install -m 644 src/mod_files/*.mod ${pkg_install_dir}/include/cuda >> install.log 2>&1
        SIRIUS_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
        cd ..
      fi
      SIRIUS_CFLAGS="-I'${pkg_install_dir}/include/cuda'"
      SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage8/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    require_env FFTW_LDFLAGS
    require_env FFTW_LIBS
    require_env FFTW_CFLAGS
    require_env ELPA_ROOT
    require_env ELPA_LDFLAGS
    require_env ELPA_LIBS
    require_env ELPA_CFLAGS
    require_env GSL_LDFLAGS
    require_env GSL_CFLAGS
    require_env GSL_LIBS
    require_env MATH_LIBS
    require_env MPI_LDFLAGS
    require_env MPI_LIBS
    require_env SCALAPACK_ROOT
    require_env SCALAPACK_LDFLAGS
    require_env SCALAPACK_CFLAGS
    require_env SCALAPACK_LIBS
    require_env LIBXC_LIBS
    require_env LIBXC_CFLAGS
    require_env LIBXC_LDFLAGS
    require_env SPGLIB_LIBS
    require_env SPGLIB_CFLAGS
    require_env SPGLIB_LDFLAGS
    require_env HDF5_LIBS
    require_env HDF5_CFLAGS
    require_env HDF5_LDFLAGS
    require_env LIBVDWXC_CFLAGS
    require_env LIBVDWXC_LDFLAGS
    require_env LIBVDWXC_LIBS
    require_env SPFFT_ROOT
    require_env SPFFT_CFLAGS
    require_env SPFFT_LDFLAGS
    require_env SPFFT_LIBS
    require_env SPLA_ROOT
    require_env SPLA_CFLAGS
    require_env SPLA_LDFLAGS
    require_env SPLA_LIBS
    check_lib -lsirius "sirius"
    add_include_from_paths SIRIUS_CFLAGS "sirius*" $INCLUDE_PATHS
    add_lib_from_paths SIRIUS_LDFLAGS "libsirius.*" $LIB_PATHS
    ;;
  *)
    echo "==================== Linking SIRIUS_Dist to user paths ===================="
    pkg_install_dir="$with_sirius"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/lib64"
    check_dir "${pkg_install_dir}/include"
    ;;
esac
if [ "$with_sirius" != "__DONTUSE__" ]; then
  SIRIUS_LIBS="-lsirius IF_CUDA(-lcusolver|)"
  SIRIUS_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
  SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
  SIRIUS_CFLAGS="-I'${pkg_install_dir}/include'"
  if [ "$with_sirius" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_sirius"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib/cuda"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib/cuda"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib/cuda"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_sirius" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_sirius"
export SIRIUS_CFLAGS="IF_CUDA(-I${pkg_install_dir}/include/cuda|-I${pkg_install_dir}/include)"
export SIRIUS_FFLAGS="IF_CUDA(-I${pkg_install_dir}/include/cuda|-I${pkg_install_dir}/include)"
export SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
export SIRIUS_CUDA_LDFLAGS="-L'${pkg_install_dir}/lib/cuda' -Wl,-rpath='${pkg_install_dir}/lib/cuda'"
export SIRIUS_LIBS="${SIRIUS_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI("-D__SIRIUS"|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI("\${SIRIUS_CFLAGS}"|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(IF_CUDA("\${SIRIUS_CUDA_LDFLAGS}"|"\${SIRIUS_LDFLAGS}")|)"
export CP_LIBS="IF_MPI("\${SIRIUS_LIBS}"|) \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_sirius"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "sirius"
