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
    add_include_from_paths CUSOLVERMP_CFLAGS "cusolverMp.h" $INCLUDE_PATHS
    add_lib_from_paths CUSOLVERMP_LDFLAGS "libcusolverMp.*" $LIB_PATHS
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
    CUSOLVERMP_CFLAGS="-I'${pkg_install_dir}/include'"
    CUSOLVERMP_LDFLAGS="-L'${CUSOLVERMP_LIBDIR}' -Wl,-rpath,'${CUSOLVERMP_LIBDIR}'"
    cusolvermp_header="${pkg_install_dir}/include/cusolverMp.h"
    cusolvermp_lib="$(find_in_paths "libcusolverMp.*" CUSOLVERMP_LIBDIR)"
    ;;
esac

if [ "${with_cusolvermp}" != "__DONTUSE__" ]; then
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

  CUSOLVERMP_LIBS="-lcusolverMp"
  CUSOLVERMP_DFLAGS="-D__CUSOLVERMP"
  if [ "${cusolvermp_major}" -gt 0 ] || [ "${cusolvermp_minor}" -ge 7 ]; then
    add_include_from_paths NCCL_CFLAGS "nccl.h" $INCLUDE_PATHS
    add_lib_from_paths NCCL_LDFLAGS "libnccl.*" $LIB_PATHS
    nccl_lib="$(find_in_paths "libnccl.*" $LIB_PATHS)"
    if [ "${nccl_lib}" = "__FALSE__" ]; then
      report_error ${LINENO} "Could not find NCCL required by cuSOLVERMp ${cusolvermp_major}.${cusolvermp_minor}."
      exit 1
    fi
    NCCL_ROOT="$(dirname "$(dirname "${nccl_lib}")")"
    CUSOLVERMP_CFLAGS="${CUSOLVERMP_CFLAGS} ${NCCL_CFLAGS}"
    CUSOLVERMP_LIBS="${CUSOLVERMP_LIBS} -lnccl"
    CUSOLVERMP_LDFLAGS="${CUSOLVERMP_LDFLAGS} ${NCCL_LDFLAGS}"
    CUSOLVERMP_DFLAGS="${CUSOLVERMP_DFLAGS} -D__CUSOLVERMP_NCCL"
  else
    add_include_from_paths CAL_CFLAGS "cal.h" $INCLUDE_PATHS
    add_lib_from_paths CAL_LDFLAGS "libcal.*" $LIB_PATHS
    cal_lib="$(find_in_paths "libcal.*" $LIB_PATHS)"
    if [ "${cal_lib}" = "__FALSE__" ]; then
      report_error ${LINENO} "Could not find CAL required by cuSOLVERMp ${cusolvermp_major}.${cusolvermp_minor}."
      exit 1
    fi
    CAL_ROOT="$(dirname "$(dirname "${cal_lib}")")"
    CUSOLVERMP_CFLAGS="${CUSOLVERMP_CFLAGS} ${CAL_CFLAGS}"
    CUSOLVERMP_LIBS="${CUSOLVERMP_LIBS} -lcal"
    CUSOLVERMP_LDFLAGS="${CUSOLVERMP_LDFLAGS} ${CAL_LDFLAGS}"
  fi

  add_include_from_paths UCC_CFLAGS "ucc/api/ucc.h" $INCLUDE_PATHS
  add_lib_from_paths UCC_LDFLAGS "libucc.*" $LIB_PATHS
  add_lib_from_paths UCX_LDFLAGS "libucs.*" $LIB_PATHS
  ucc_lib="$(find_in_paths "libucc.*" $LIB_PATHS)"
  ucx_lib="$(find_in_paths "libucs.*" $LIB_PATHS)"
  if [ "${ucc_lib}" = "__FALSE__" ] || [ "${ucx_lib}" = "__FALSE__" ]; then
    report_error ${LINENO} "Could not find UCC/UCX required by cuSOLVERMp."
    exit 1
  fi
  UCC_ROOT="$(dirname "$(dirname "${ucc_lib}")")"
  UCX_ROOT="$(dirname "$(dirname "${ucx_lib}")")"
  CUSOLVERMP_LDFLAGS="${CUSOLVERMP_LDFLAGS} ${UCC_LDFLAGS} ${UCX_LDFLAGS}"
  CUSOLVERMP_LIBS="${CUSOLVERMP_LIBS} -lucc -lucs"

  cat << EOF > "${BUILDDIR}/setup_cusolvermp"
export CUSOLVERMP_VER="${cusolvermp_major}.${cusolvermp_minor}"
export CUSOLVERMP_CFLAGS="${CUSOLVERMP_CFLAGS}"
export CUSOLVERMP_LDFLAGS="${CUSOLVERMP_LDFLAGS}"
export CUSOLVERMP_LIBS="${CUSOLVERMP_LIBS}"
export CUSOLVER_MP_ROOT="${pkg_install_dir}"
export UCC_ROOT="${UCC_ROOT}"
export UCX_ROOT="${UCX_ROOT}"
export CP_DFLAGS="\${CP_DFLAGS} IF_CUDA(${CUSOLVERMP_DFLAGS}|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_CUDA(${CUSOLVERMP_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_CUDA(${CUSOLVERMP_LDFLAGS}|)"
export CP_LIBS="IF_CUDA(${CUSOLVERMP_LIBS}|) \${CP_LIBS}"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
prepend_path CMAKE_PREFIX_PATH "${UCC_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${UCX_ROOT}"
EOF
  if [ -n "${NCCL_ROOT}" ]; then
    cat << EOF >> "${BUILDDIR}/setup_cusolvermp"
export NCCL_ROOT="${NCCL_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${NCCL_ROOT}"
EOF
  fi
  if [ -n "${CAL_ROOT}" ]; then
    cat << EOF >> "${BUILDDIR}/setup_cusolvermp"
export CAL_ROOT="${CAL_ROOT}"
prepend_path CMAKE_PREFIX_PATH "${CAL_ROOT}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_cusolvermp" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_cusolvermp"
write_toolchain_env "${INSTALLDIR}"

report_timing "cusolvermp"

#EOF
