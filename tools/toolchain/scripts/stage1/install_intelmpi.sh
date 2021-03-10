#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "intelmpi" ] && exit 0
rm -f "${BUILDDIR}/setup_intelmpi"

INTELMPI_CFLAGS=''
INTELMPI_LDFLAGS=''
INTELMPI_LIBS=''
mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_intelmpi" in
  __INSTALL__)
    echo '__INSTALL__ is not supported; please manually install Intel MPI'
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding Intel MPI from system paths ===================="
    check_command mpirun "intelmpi"
    check_command mpiicc "intelmpi"
    check_command mpiifort "intelmpi"
    check_command mpiicpc "intelmpi"
    add_include_from_paths INTELMPI_CFLAGS "mpi.h" $INCLUDE_PATHS
    add_lib_from_paths INTELMPI_LDFLAGS "libmpi.*" $LIB_PATHS
    check_lib -lmpi "intelmpi"
    check_lib -lmpicxx "intelmpi"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking INTELMPI to user paths ===================="
    pkg_install_dir="$with_intelmpi"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    INTELMPI_CFLAGS="-I'${pkg_install_dir}/include'"
    INTELMPI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_intelmpi" != "__DONTUSE__" ]; then
  INTELMPI_LIBS="-lmpi -lmpicxx"
  if [ "$with_intelmpi" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_intelmpi"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_intelmpi" >> $SETUPFILE
    mpi_bin="$pkg_install_dir/bin/mpirun"
  else
    mpi_bin=mpirun
  fi
  cat << EOF >> "${BUILDDIR}/setup_intelmpi"
export MPI_MODE="${MPI_MODE}"
export INTELMPI_CFLAGS="${INTELMPI_CFLAGS}"
export INTELMPI_LDFLAGS="${INTELMPI_LDFLAGS}"
export INTELMPI_LIBS="${INTELMPI_LIBS}"
export MPI_CFLAGS="${INTELMPI_CFLAGS}"
export MPI_LDFLAGS="${INTELMPI_LDFLAGS}"
export MPI_LIBS="${INTELMPI_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel ${mpi2_dflags}|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${INTELMPI_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${INTELMPI_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_MPI(${INTELMPI_LIBS}|)"
EOF
fi

load "${BUILDDIR}/setup_intelmpi"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "intelmpi"
