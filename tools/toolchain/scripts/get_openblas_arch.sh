#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

openblas_ver="0.3.28" # Keep in sync with install_openblas.sh
openblas_sha256="f1003466ad074e9b0c8d421a204121100b0751c96fc6fcf3d1456bd12f8a00a1"
openblas_pkg="OpenBLAS-${openblas_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

find_openblas_dir() {
  local __dir=''
  for __dir in *OpenBLAS*; do
    if [ -d "$__dir" ]; then
      echo "$__dir"
      return 0
    fi
  done
  echo ''
}

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

echo "==================== Getting proc arch info using OpenBLAS tools ===================="
# find existing openblas source dir
openblas_dir="$(find_openblas_dir)"
# if cannot find openblas source dir, try download one
if ! [ "$openblas_dir" ]; then
  if [ -f ${openblas_pkg} ]; then
    echo "${openblas_pkg} is found"
  else
    download_pkg_from_cp2k_org "${openblas_sha256}" "${openblas_pkg}"
  fi
  tar -xzf ${openblas_pkg}
  openblas_dir="$(find_openblas_dir)"
fi
openblas_conf="${openblas_dir}/Makefile.conf"
# try find Makefile.config, if not then generate one with make lapack_prebuild
if ! [ -f "$openblas_conf" ]; then
  cd "$openblas_dir"
  make lapack_prebuild
  cd ..
fi
OPENBLAS_LIBCORE="$(grep 'LIBCORE=' $openblas_conf | cut -f2 -d=)"
OPENBLAS_ARCH="$(grep 'ARCH=' $openblas_conf | cut -f2 -d=)"
echo "OpenBLAS detected LIBCORE = $OPENBLAS_LIBCORE"
echo "OpenBLAS detected ARCH    = $OPENBLAS_ARCH"
# output setup file
cat << EOF > "${BUILDDIR}/openblas_arch"
export OPENBLAS_LIBCORE="${OPENBLAS_LIBCORE}"
export OPENBLAS_ARCH="${OPENBLAS_ARCH}"
EOF
