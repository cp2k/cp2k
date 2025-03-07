#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

dbcsr_ver="2.8.0"
dbcsr_sha256="d55e4f052f28d1ed0faeaa07557241439243287a184d1fd27f875c8b9ca6bd96"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_dbcsr" ] && rm "${BUILDDIR}/setup_dbcsr"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_dbcsr}" in
  __INSTALL__)
    echo "==================== Installing DBCSR ===================="
    pkg_install_dir="${INSTALLDIR}/dbcsr-${dbcsr_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "dbcsr-${dbcsr_ver} is already installed, skipping it."
    else
      if [ -f dbcsr-${dbcsr_ver}.tar.gz ]; then
        echo "dbcsr-${dbcsr_ver}.tar.gz is found"
      else
        download_pkg_from_urlpath "${dbcsr_sha256}" "dbcsr-${dbcsr_ver}.tar.gz" \
          https://github.com/cp2k/dbcsr/releases/download/v${dbcsr_ver}
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d dbcsr-${dbcsr_ver} ] && rm -rf dbcsr-${dbcsr_ver}
      tar -xzf dbcsr-${dbcsr_ver}.tar.gz
      cd dbcsr-${dbcsr_ver}
      mkdir build-cpu
      cd build-cpu
      CMAKE_OPTIONS="-DBUILD_TESTING=NO -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_VERBOSE_MAKEFILE=ON"
      CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DUSE_OPENMP=ON -DUSE_SMM=blas -DWITH_EXAMPLES=NO"
      if [ "${MPI_MODE}" == "no" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_MPI=OFF"
      else
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_MPI=ON"
      fi
      cmake \
        -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        ${CMAKE_OPTIONS} .. \
        > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      if [ "${ENABLE_CUDA}" == "__TRUE__" ]; then
        mkdir build-cuda
        cd build-cuda
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_ACCEL=cuda -DWITH_GPU=P100"
        cmake \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir}-cuda \
          ${CMAKE_OPTIONS} .. \
          > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi
      if [ "${ENABLE_HIP}" == "__TRUE__" ]; then
        mkdir build-hip
        cd build-hip
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_ACCEL=hip -DWITH_GPU=Mi250"
        cmake \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir}-hip \
          ${CMAKE_OPTIONS} .. \
          > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
        make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
        make -j $(get_nprocs) install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
        cd ..
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage9/$(basename ${SCRIPT_NAME})"
      DBCSR_CFLAGS="-I'${pkg_install_dir}/include'"
      DBCSR_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
      DBCSR_CUDA_CFLAGS="-I'${pkg_install_dir}-cuda/include'"
      DBCSR_CUDA_LDFLAGS="-L'${pkg_install_dir}-cuda/lib' -Wl,-rpath='${pkg_install_dir}-cuda/lib'"
      DBCSR_HIP_CFLAGS="-I'${pkg_install_dir}-hip/include'"
      DBCSR_HIP_LDFLAGS="-L'${pkg_install_dir}-hip/lib' -Wl,-rpath='${pkg_install_dir}-hip/lib'"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding DBCSR from system paths ===================="
    check_lib -ldbcsr "dbcsr"
    add_include_from_paths DBCSR_CFLAGS "dbcsr.h" $INCLUDE_PATHS
    add_lib_from_paths DBCSR_LDFLAGS "dbcsr.*" $LIB_PATHS
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking DBCSR to user paths ===================="
    pkg_install_dir="${with_dbcsr}"
    DBCSR_LIBDIR="${pkg_install_dir}/lib"
    check_dir "${DBCSR_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    DBCSR_CFLAGS="-I'${pkg_install_dir}/include'"
    DBCSR_LDFLAGS="-L'${DBCSR_LIBDIR}' -Wl,-rpath,'${DBCSR_LIBDIR}'"
    ;;
esac

if [ "${with_dbcsr}" != "__DONTUSE__" ]; then
  DBCSR_LIBS="-ldbcsr"
  if [ "${with_dbcsr}" != "__SYSTEM__" ]; then
    if [ "${ENABLE_CUDA}" == "__TRUE__" ]; then
      pkg_install_dir1="${pkg_install_dir}-cuda"
    else
      if [ "${ENABLE_HIP}" == "__TRUE__" ]; then
        pkg_install_dir1="${pkg_install_dir}-hip"
      else
        pkg_install_dir1="${pkg_install_dir}"
      fi
    fi
  fi
  cat << EOF > "${BUILDDIR}/setup_dbcsr"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir1}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir1}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir1}/lib"
prepend_path CPATH "${pkg_install_dir1}/include"
prepend_path CMAKE_INSTALL_PREFIX "${pkg_install_dir1}"
export DBCSR_ROOT="${pkg_install_dir}"
export DBCSR_HIP_ROOT="${pkg_install_dir}-hip"
export DBCSR_CUDA_ROOT="${pkg_install_dir}-cuda"
export DBCSR_VER="${dbcsr_ver}"
export DBCSR_DIR="${pkg_install_dir1}/lib/cmake/dbcsr"
export DBCSR_CFLAGS="${DBCSR_CFLAGS}"
export DBCSR_LDFLAGS="IF_CUDA(${DBCSR_CUDA_LDFLAGS}|IF_HIP(${DBCSR_HIP_LDFLAGS}|${DBCSR_LDFLAGS}))"
export DBCSR_LIBS="${DBCSR_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_CUDA(-D__DBCSR_ACC -D__DBCSR|IF_HIP(-D__DBCSR_ACC -D__DBCSR|-D__DBCSR))"
export CP_CFLAGS="\${CP_CFLAGS} ${DBCSR_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${DBCSR_LDFLAGS}"
export CP_LIBS="${DBCSR_LIBS} \${CP_LIBS}"
EOF
else
  touch "${BUILDDIR}/setup_dbcsr"
fi

cat "${BUILDDIR}/setup_dbcsr" >> ${SETUPFILE}

load "${BUILDDIR}/setup_dbcsr"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "DBCSR"
