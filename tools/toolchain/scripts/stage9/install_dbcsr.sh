#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

dbcsr_ver="2.10.0"
dbcsr_sha256="3d897220fbb4498215331efad6905eb7744881b4cf04eb5c5fb4db7c48a56ef9"
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
      # retrieve_package "${dbcsr_sha256}" "dbcsr-${dbcsr_ver}.tar.gz"
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
      # DBCSR may predate GB10. Build native sm_121 code while
      # reusing the closest available libsmm_acc parameters.
      if [ "${ENABLE_CUDA}" == "__TRUE__" ] && [ "${GPUVER}" == "GB10" ]; then
        if ! grep -q "GB10" CMakeLists.txt; then
          sed -i "s/    H100)/    H100\\n    GB10)/" CMakeLists.txt
          sed -i "/  set(GPU_ARCH_NUMBER_H100 90)/a\\  set(GPU_ARCH_NUMBER_GB10 121)" CMakeLists.txt
          cp src/acc/libsmm_acc/parameters/parameters_H100.json \
            src/acc/libsmm_acc/parameters/parameters_GB10.json
        fi
      fi
      mkdir build-cpu
      cd build-cpu
      CMAKE_OPTIONS="-DBUILD_TESTING=NO -DCMAKE_INSTALL_LIBDIR=lib -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_VERBOSE_MAKEFILE=ON"
      CMAKE_OPTIONS="${CMAKE_OPTIONS} -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DUSE_OPENMP=ON -DWITH_EXAMPLES=NO"
      if [ "${with_libxs}" != "__DONTUSE__" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_LIBXS=ON"
        if [ "${with_libxsmm}" != "__DONTUSE__" ]; then
          CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_LIBXSMM=ON"
        fi
      fi
      if [ "${MPI_MODE}" == "no" ]; then
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_MPI=OFF"
      else
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_MPI=ON -DUSE_MPI_F08=ON"
      fi
      cmake \
        -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        ${CMAKE_OPTIONS} .. \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
      cd ..
      if [ "${ENABLE_CUDA}" == "__TRUE__" ]; then
        echo "Installing from scratch into ${pkg_install_dir}-cuda"
        mkdir build-cuda
        cd build-cuda
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_ACCEL=cuda"
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DWITH_GPU=${GPUVER:-P100}"
        cmake \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir}-cuda \
          ${CMAKE_OPTIONS} .. \
          > cmake.log 2>&1 || tail_excerpt cmake.log
        make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
        cd ..
      fi
      if [ "${ENABLE_HIP}" == "__TRUE__" ]; then
        echo "Installing from scratch into ${pkg_install_dir}-hip"
        mkdir build-hip
        cd build-hip
        CMAKE_OPTIONS="${CMAKE_OPTIONS} -DUSE_ACCEL=hip -DWITH_GPU=Mi250"
        cmake \
          -DCMAKE_INSTALL_PREFIX=${pkg_install_dir}-hip \
          ${CMAKE_OPTIONS} .. \
          > cmake.log 2>&1 || tail_excerpt cmake.log
        make -j $(get_nprocs) install > make.log 2>&1 || tail_excerpt make.log
        cd ..
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage9/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding DBCSR from system paths ===================="
    check_lib -ldbcsr "dbcsr"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libdbcsr.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking DBCSR to user paths ===================="
    pkg_install_dir="${with_dbcsr}"
    DBCSR_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && DBCSR_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${DBCSR_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac

if [ "${with_dbcsr}" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_dbcsr"
export DBCSR_VER="${dbcsr_ver}"
EOF
  if [ "${with_dbcsr}" != "__SYSTEM__" ]; then
    pkg_install_dir1="${pkg_install_dir}"
    if [ "${with_dbcsr}" = "__INSTALL__" ]; then
      if [ "${ENABLE_CUDA}" = "__TRUE__" ]; then
        pkg_install_dir1="${pkg_install_dir}-cuda"
      elif [ "${ENABLE_HIP}" = "__TRUE__" ]; then
        pkg_install_dir1="${pkg_install_dir}-hip"
      fi
    fi
    cat << EOF >> "${BUILDDIR}/setup_dbcsr"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir1}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir1}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir1}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir1}"
EOF
  fi
else
  touch "${BUILDDIR}/setup_dbcsr"
fi

filter_setup "${BUILDDIR}/setup_dbcsr" "${SETUPFILE}"

load "${BUILDDIR}/setup_dbcsr"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "DBCSR"
