#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=${0}
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_amd" ] && rm "${BUILDDIR}/setup_amd"

mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_amd}" in
  __INSTALL__)
    echo "==================== Installing the AMD compiler ======================"
    echo "Installation of AMD compilers is not supported; please install manually"
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding AMD compiler from system paths ===================="
    check_command clang "amd" && CC="$(real_path $(command -v clang))" || exit 1
    check_command clang++ "amd" && CXX="$(real_path $(command -v clang++))" || exit 1
    check_command flang "amd" && FC="$(real_path $(command -v flang))" || exit 1
    F90="${FC}"
    F77="${FC}"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking AMD compiler to user paths ===================="
    pkg_install_dir="${with_amd}"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_command ${pkg_install_dir}/bin/clang "amd" && CC="${pkg_install_dir}/bin/clang" || exit 1
    check_command ${pkg_install_dir}/bin/clang++ "amd" && CXX="${pkg_install_dir}/bin/clang++" || exit 1
    check_command ${pkg_install_dir}/bin/flang "amd" && FC="${pkg_install_dir}/bin/flang" || exit 1
    F90="${FC}"
    F77="${FC}"
    ;;
esac
if [ "${with_amd}" != "__DONTUSE__" ]; then
  echo "CC  is ${CC}"
  [ $(real_path $(command -v clang) | grep -e aocc-compiler) ] || echo "Check the AMD C compiler path"
  echo "CXX is ${CXX}"
  [ $(real_path $(command -v clang++) | grep -e aocc-compiler) ] || echo "Check the AMD C++ compiler path"
  echo "FC  is ${FC}"
  [ $(real_path $(command -v flang) | grep -e aocc-compiler) ] || echo "Check the AMD Fortran compiler path"
  cat << EOF > "${BUILDDIR}/setup_amd"
export CC="${CC}"
export CXX="${CXX}"
export FC="${FC}"
export F90="${F90}"
export F77="${F77}"
EOF
  if [ "${with_amd}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_amd"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_amd" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_amd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "amd"
