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

[ -f "${BUILDDIR}/setup_cusolvermp" ] && rm "${BUILDDIR}/setup_cusolvermp"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_cusolvermp}" in
  __INSTALL__)
    report_error ${LINENO} "Installing cuSOLVERMp is not supported. Use --with-cusolvermp=system or --with-cusolvermp=<prefix>."
    exit 1
    ;;
  __SYSTEM__)
    echo "==================== Finding cuSOLVERMp from system paths ===================="
    cusolvermp_header="$(find_in_paths "cusolverMp.h" $INCLUDE_PATHS)"
    cusolvermp_lib="$(find_in_paths "libcusolverMp.*" $LIB_PATHS)"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking cuSOLVERMp to user paths ===================="
    pkg_install_dir="${with_cusolvermp}"
    CUSOLVERMP_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && CUSOLVERMP_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${CUSOLVERMP_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    cusolvermp_header="${pkg_install_dir}/include/cusolverMp.h"
    cusolvermp_lib="$(find_in_paths "libcusolverMp.*" CUSOLVERMP_LIBDIR)"
    ;;
esac

if [ "${with_cusolvermp}" != "__DONTUSE__" ]; then
  # These values can be inherited from a previous toolchain run. Only export
  # roots for the communication backend selected by the detected version.
  unset NCCL_ROOT CAL_ROOT UCC_ROOT UCX_ROOT

  if [ "${cusolvermp_header}" = "__FALSE__" ] || ! [ -f "${cusolvermp_header}" ]; then
    report_error ${LINENO} "Could not find cusolverMp.h."
    exit 1
  fi
  if [ "${cusolvermp_lib}" = "__FALSE__" ] || ! [ -f "${cusolvermp_lib}" ]; then
    report_error ${LINENO} "Could not find libcusolverMp."
    exit 1
  fi

  pkg_install_dir="$(dirname "$(dirname "${cusolvermp_lib}")")"
  cusolvermp_major="$(grep -E '^#define CUSOLVERMP_VER_MAJOR' "${cusolvermp_header}" | awk '{print $3}')"
  cusolvermp_minor="$(grep -E '^#define CUSOLVERMP_VER_MINOR' "${cusolvermp_header}" | awk '{print $3}')"
  cusolvermp_major="${cusolvermp_major:-0}"
  cusolvermp_minor="${cusolvermp_minor:-0}"

  if [ "${cusolvermp_major}" -gt 0 ] || [ "${cusolvermp_minor}" -ge 7 ]; then
    cusolvermp_comm_backend="nccl"
    nccl_lib="$(find_in_paths "libnccl.*" $LIB_PATHS)"
    if [ "${nccl_lib}" = "__FALSE__" ]; then
      report_error ${LINENO} "Could not find NCCL required by cuSOLVERMp ${cusolvermp_major}.${cusolvermp_minor}."
      exit 1
    fi
    NCCL_ROOT="$(dirname "$(dirname "${nccl_lib}")")"
  else
    cusolvermp_comm_backend="cal"
    cal_lib="$(find_in_paths "libcal.*" $LIB_PATHS)"
    if [ "${cal_lib}" = "__FALSE__" ]; then
      report_error ${LINENO} "Could not find CAL required by cuSOLVERMp ${cusolvermp_major}.${cusolvermp_minor}."
      exit 1
    fi
    CAL_ROOT="$(dirname "$(dirname "${cal_lib}")")"
    ucc_lib="$(find_in_paths "libucc.*" $LIB_PATHS)"
    ucx_lib="$(find_in_paths "libucs.*" $LIB_PATHS)"
    if [ "${ucc_lib}" = "__FALSE__" ] || [ "${ucx_lib}" = "__FALSE__" ]; then
      report_error ${LINENO} "Could not find UCC/UCX required by cuSOLVERMp ${cusolvermp_major}.${cusolvermp_minor}."
      exit 1
    fi
    UCC_ROOT="$(dirname "$(dirname "${ucc_lib}")")"
    UCX_ROOT="$(dirname "$(dirname "${ucx_lib}")")"
  fi

  cat << EOF > "${BUILDDIR}/setup_cusolvermp"
export CUSOLVERMP_VER="${cusolvermp_major}.${cusolvermp_minor}"
export CUSOLVERMP_LIBS="${CUSOLVERMP_LIBS}"
export CUSOLVER_MP_ROOT="${pkg_install_dir}"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  if [ "${cusolvermp_comm_backend}" = "nccl" ]; then
    cat << EOF >> "${BUILDDIR}/setup_cusolvermp"
export NCCL_ROOT="${NCCL_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${NCCL_ROOT}"
EOF
  else
    cat << EOF >> "${BUILDDIR}/setup_cusolvermp"
export CAL_ROOT="${CAL_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${CAL_ROOT}"
export UCC_ROOT="${UCC_ROOT}"
export UCX_ROOT="${UCX_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${UCC_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${UCX_ROOT}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_cusolvermp" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_cusolvermp"
write_toolchain_env "${INSTALLDIR}"

report_timing "cusolvermp"

#EOF
