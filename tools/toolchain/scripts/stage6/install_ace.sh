#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

ace_ver="v.2023.11.25.fix2"
ace_dir="lammps-user-pace-${ace_ver}"
ace_pkg="${ace_dir}.tar.gz"
ace_sha256="e0885351a8a730f5576dace2374fa470523a4526383c6a64af571e1344a40686"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_ace" ] && rm "${BUILDDIR}/setup_ace"

ACE_LDFLAGS=''
ACE_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_ace" in
  __INSTALL__)
    echo "==================== Installing Ace ======================="
    pkg_install_dir="${INSTALLDIR}/${ace_dir}"
    install_lock_file="${pkg_install_dir}/install_successful"
    ace_root="${pkg_install_dir}"
    if verify_checksums "${install_lock_file}"; then
      echo "${ace_dir} aka Ace is already installed, skipping it."
    else
      if [ -f ${ace_pkg} ]; then
        echo "${ace_pkg} is found"
      else
        download_pkg_from_urlpath "${ace_sha256}" "${ace_ver}.tar.gz" "https://github.com/ICAMS/lammps-user-pace/archive/refs/tags" "${ace_pkg}"
      fi
      [ -d ${ace_dir} ] && rm -rf ${ace_dir}
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf ${ace_pkg}
      cd ${ace_dir}

      # Fix for GCC 15
      sed -i '1i #include <cstdint>' yaml-cpp/src/emitterutils.cpp

      mkdir build
      cd build

      cmake \
        -DCMAKE_CXX_STANDARD=17 \
        .. > cmake.log 2>&1 || tail -n ${LOG_LINES} cmake.log
      make -j > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      # no make install.
      [ -d ${pkg_install_dir} ] && rm -rf ${pkg_install_dir}
      mkdir -p ${pkg_install_dir}/lib
      cp -a libpace.a libcnpy.a build-yaml-cpp/libyaml-cpp-pace.a \
        ${pkg_install_dir}/lib
      cp -a ../yaml-cpp/include ${pkg_install_dir}
      mkdir ${pkg_install_dir}/include/ace
      cp -a ../ML-PACE/ace/*.h ${pkg_install_dir}/include/ace
      mkdir ${pkg_install_dir}/include/ace-evaluator
      cp -a ../ML-PACE/ace-evaluator/*.h ${pkg_install_dir}/include/ace-evaluator
      #
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    ACE_CFLAGS="-I'${pkg_install_dir}/include'"
    #smuggle include dirs to CXXFLAGS via DFLAGS....
    ACE_DFLAGS="-D__ACE ${ACE_CFLAGS}"
    ACE_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
    #    not supported
    #  __SYSTEM__)
    #    echo "==================== Finding Ace from system paths ===================="
    #    check_lib -lace "ACE"
    #    add_lib_from_paths ACE_LDFLAGS "libpace*" $LIB_PATHS
    #    add_include_from_paths ACE_CFLAGS "ace" $INCLUDE_PATHS
    #    ACE_DFLAGS="-D__ACE"
    #    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking ACE to user paths ===================="
    pkg_install_dir="$with_ace"
    check_dir "${pkg_install_dir}/include/ace"
    check_dir "${pkg_install_dir}/include/ace-evaluator"
    check_dir "${pkg_install_dir}/include/yaml-cpp"
    check_dir "${pkg_install_dir}/lib"
    ACE_CFLAGS="-I'${pkg_install_dir}/include'"
    #smuggle include dirs to CXXFLAGS via DFLAGS....
    ACE_DFLAGS="-D__ACE ${ACE_CFLAGS}"
    ACE_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_ace" != "__DONTUSE__" ]; then
  ACE_LIBS='-Wl,--start-group -lpace -lyaml-cpp-pace -lcnpy -Wl,--end-group -lstdc++'
  cat << EOF > "${BUILDDIR}/setup_ace"
export ACE_VER="${ace_ver}"
EOF
  if [ "$with_ace" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_ace"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
    cat "${BUILDDIR}/setup_ace" >> $SETUPFILE
  fi

  cat << EOF >> "${BUILDDIR}/setup_ace"
export ACE_DFLAGS="${ACE_DFLAGS}"
export ACE_CFLAGS="${ACE_CFLAGS}"
export ACE_LDFLAGS="${ACE_LDFLAGS}"
export ACE_LIBS="${ACE_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} ${ACE_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${ACE_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${ACE_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${ACE_LIBS}"
EOF
#  cat << EOF >> "${INSTALLDIR}/lsan.supp"
## leaks related to ACE
#EOF
fi

load "${BUILDDIR}/setup_ace"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "ace"
