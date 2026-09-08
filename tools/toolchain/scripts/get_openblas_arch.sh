#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

openblas_ver="0.3.34" # Keep in sync with install_openblas.sh
openblas_sha256="cd7e129868320cc2d033afa920e31202dfe0b8066a5b66661900ccc0f197dfed"
openblas_pkg="OpenBLAS-${openblas_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

find_openblas_dir() {
  find . -maxdepth 1 -type d -name '*OpenBLAS*' 2> /dev/null | head -1
}

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

echo "==================== Getting proc arch info using OpenBLAS tools ===================="

get_system_openblas_arch() {
  local openblas_lib=""
  local openblas_dir=""

  if [ "${with_openblas}" = "__SYSTEM__" ]; then
    openblas_lib=$(find /usr/lib /usr/local/lib -maxdepth 3 -name "libopenblas*" -type f 2> /dev/null | head -1)
  elif [ "${with_openblas}" != "__INSTALL__" ] && [ -d "${with_openblas}" ]; then
    openblas_lib=$(find "${with_openblas}/lib" -maxdepth 2 -name "libopenblas*" -type f 2> /dev/null | head -1)
  fi

  if [ -z "$openblas_lib" ]; then
    return 1
  fi

  openblas_dir="$(dirname "$(dirname "$openblas_lib")")"

  local makefile_conf=""
  for dir in "${openblas_dir}/include" "${openblas_dir}/share/OpenBLAS" "/usr/share/OpenBLAS" "/usr/include"; do
    if [ -f "${dir}/Makefile.conf" ]; then
      makefile_conf="${dir}/Makefile.conf"
      break
    fi
  done

  if [ -n "$makefile_conf" ]; then
    OPENBLAS_LIBCORE="$(sed -n 's/^LIBCORE=//p' "$makefile_conf")"
    OPENBLAS_ARCH="$(sed -n 's/^ARCH=//p' "$makefile_conf")"
    return 0
  fi

  local openblas_so="$(readelf -A "$openblas_lib" 2> /dev/null | grep 'Unknown' | head -1)"
  case "$openblas_so" in
    *x86_64*) OPENBLAS_ARCH="x86_64" ;;
    *aarch64* | *arm64*) OPENBLAS_ARCH="arm64" ;;
    *ppc64*) OPENBLAS_ARCH="ppc64" ;;
    *) OPENBLAS_ARCH="$(uname -m)" ;;
  esac
  OPENBLAS_LIBCORE=""
  return 0
}

if get_system_openblas_arch; then
  echo "Using system OpenBLAS architecture detection"
else
  openblas_dir="$(find_openblas_dir)"
  if [ -z "$openblas_dir" ]; then
    retrieve_package "${openblas_sha256}" "${openblas_pkg}"
    tar -xzf "${openblas_pkg}"
    openblas_dir="$(find_openblas_dir)"
  fi
  openblas_conf="${openblas_dir}/Makefile.conf"
  if ! [ -f "$openblas_conf" ]; then
    cd "$openblas_dir"
    make lapack_prebuild
    cd ..
  fi
  OPENBLAS_LIBCORE="$(sed -n 's/^LIBCORE=//p' "$openblas_conf")"
  OPENBLAS_ARCH="$(sed -n 's/^ARCH=//p' "$openblas_conf")"
fi

echo "OpenBLAS detected LIBCORE = $OPENBLAS_LIBCORE"
echo "OpenBLAS detected ARCH    = $OPENBLAS_ARCH"
cat << EOF > "${BUILDDIR}/openblas_arch"
export OPENBLAS_LIBCORE="${OPENBLAS_LIBCORE}"
export OPENBLAS_ARCH="${OPENBLAS_ARCH}"
EOF
