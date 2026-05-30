#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

ace_ver="2025.12.4.p1"
ace_dir="lammps-user-pace-v.${ace_ver}"
ace_pkg="ace-${ace_ver}.tar.gz"
ace_sha256="21e9d7ad2094eef0f19958d154866fc725fc6ccfa82ec3681ef2b006545ced96"

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
      retrieve_package "${ace_sha256}" "${ace_pkg}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d ${ace_dir} ] && rm -rf ${ace_dir}
      tar -xzf ${ace_pkg}
      cd ${ace_dir}

      # Fix for GCC 15 (portable sed: prepend line 1 - Works on macOS + Linux; Doesn’t depend on sed’s i\ quirks; Avoids adding the include twice if you rebuild)
      if ! grep -q '^#include <cstdint>' yaml-cpp/src/emitterutils.cpp; then
        if [ "$(uname)" = "Darwin" ]; then
          sed -i '' '1s;^;#include <cstdint>\n;' yaml-cpp/src/emitterutils.cpp
        else
          sed -i '1s;^;#include <cstdint>\n;' yaml-cpp/src/emitterutils.cpp
        fi
      fi

      rm -rf build
      mkdir -p build && cd build
      # fix: without DCMAKE_DISABLE_FIND_PACKAGE_yaml-cpp the cp line below will crash in all those cases where yaml is system installed.
      cmake .. \
        -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_DISABLE_FIND_PACKAGE_yaml-cpp=TRUE \
        > cmake.log 2>&1 || tail_excerpt cmake.log
      make install -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    ACE_CFLAGS="-I'${pkg_install_dir}/include'"
    #smuggle include dirs to CXXFLAGS via DFLAGS....
    ACE_DFLAGS="-D__ACE ${ACE_CFLAGS}"
    ACE_LDFLAGS="-L'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding Ace from system paths ===================="
    # check_lib -lace "ACE"
    # add_lib_from_paths ACE_LDFLAGS "libpace*" $LIB_PATHS
    # add_include_from_paths ACE_CFLAGS "ace" $INCLUDE_PATHS
    # ACE_DFLAGS="-D__ACE"
    report_error "The system detection of Ace is not supported; please use \"--with-ace=install\"."
    exit 1
    ;;
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
    filter_setup "${BUILDDIR}/setup_ace" "${SETUPFILE}"
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
