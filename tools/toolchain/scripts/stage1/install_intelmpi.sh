#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "intelmpi" ] && exit 0
rm -f "${BUILDDIR}/setup_intelmpi"

mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_intelmpi}" in
  __INSTALL__)
    echo "==================== Installing Intel MPI ===================="
    echo "Installation of Intel MPI is not supported; please install manually"
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding Intel MPI from system paths ===================="
    check_command mpiexec "intelmpi" && MPIEXEC="$(real_path $(command -v mpiexec))" || exit 1
    if [ "${with_intel}" != "__DONTUSE__" ]; then
      if [ "${with_ifx}" = "yes" ]; then
        check_command mpiicx "intelmpi" && MPICC="$(real_path $(command -v mpiicx))" || exit 1
        check_command mpiicpx "intelmpi" && MPICXX="$(real_path $(command -v mpiicpx))" || exit 1
        check_command mpiifx "intelmpi" && MPIFC="$(real_path $(command -v mpiifx))" || exit 1
      else
        check_command mpiicc "intelmpi" && MPICC="$(real_path $(command -v mpiicc))" || exit 1
        check_command mpiicpc "intelmpi" && MPICXX="$(real_path $(command -v mpiicpc))" || exit 1
        check_command mpiifort "intelmpi" && MPIFC="$(real_path $(command -v mpiifort))" || exit 1
      fi
    else
      echo "The use of Intel MPI is only supported with the Intel compiler"
      exit 1
    fi
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    check_lib -lmpi "intelmpi"
    check_lib -lmpicxx "intelmpi"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking Intel MPI to user paths ===================="
    pkg_install_dir="${with_intelmpi}"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_command ${pkg_install_dir}/bin/mpiexec "intel" && MPIEXEC="${pkg_install_dir}/bin/mpiexec" || exit 1
    if [ "${with_intel}" != "__DONTUSE__" ]; then
      check_command ${pkg_install_dir}/bin/mpiicc "intel" && MPICC="${pkg_install_dir}/bin/mpiicc" || exit 1
      check_command ${pkg_install_dir}/bin/mpiicpc "intel" && MPICXX="${pkg_install_dir}/bin/mpiicpc" || exit 1
      check_command ${pkg_install_dir}/bin/mpiifort "intel" && MPIFC="${pkg_install_dir}/bin/mpiifort" || exit 1
    else
      echo "The use of Intel MPI is only supported with the Intel compiler"
      exit 1
    fi
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    ;;
esac
if [ "${with_intelmpi}" != "__DONTUSE__" ]; then
  I_MPI_CXX="icpx"
  I_MPI_CC="icx"
  if [ "${with_ifx}" = "yes" ]; then
    I_MPI_FC="ifx"
  else
    I_MPI_FC="ifort"
  fi
  echo "I_MPI_CXX is ${I_MPI_CXX}"
  echo "I_MPI_CC  is ${I_MPI_CC}"
  echo "I_MPI_FC  is ${I_MPI_FC}"
  echo "MPICXX    is ${MPICXX}"
  echo "MPICC     is ${MPICC}"
  echo "MPIFC     is ${MPIFC}"
  cat << EOF > "${BUILDDIR}/setup_intelmpi"
export I_MPI_CXX="${I_MPI_CXX}"
export I_MPI_CC="${I_MPI_CC}"
export I_MPI_FC="${I_MPI_FC}"
export MPI_MODE="${MPI_MODE}"
export MPIEXEC="${MPIEXEC}"
export MPICC="${MPICC}"
export MPICXX="${MPICXX}"
export MPIFC="${MPIFC}"
export MPIFORT="${MPIFORT}"
export MPIF77="${MPIF77}"
EOF
  if [ "${with_intelmpi}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_intelmpi"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_intelmpi" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_intelmpi"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "intelmpi"
