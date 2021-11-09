#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

openblas_ver="0.3.18" # Keep in sync with install_openblas.sh
openblas_sha256="1632c1e8cca62d8bed064b37747e331a1796fc46f688626337362bf0d16aeadb"
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
    download_pkg ${DOWNLOADER_FLAGS} ${openblas_sha256} \
      https://www.cp2k.org/static/downloads/${openblas_pkg}
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
