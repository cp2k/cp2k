#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

find_openblas_dir() {
    local __dir=''
    for __dir in *OpenBLAS* ; do
        if [ -d "$__dir" ] ; then
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
if ! [ "$openblas_dir" ] ; then
    if [ -f openblas-${openblas_ver}.tar.gz ] ; then
        echo "openblas-${openblas_ver}.tar.gz is found"
    else
        download_pkg ${DOWNLOADER_FLAGS} \
                     https://www.cp2k.org/static/downloads/OpenBLAS-${openblas_ver}.tar.gz
    fi
    tar -xzf OpenBLAS-${openblas_ver}.tar.gz
    openblas_dir="$(find_openblas_dir)"
fi
openblas_conf="${openblas_dir}/Makefile.conf"
# try find Makefile.config, if not then generate one with make lapack_prebuild
if ! [ -f "$openblas_conf" ] ; then
    cd "$openblas_dir"
    make lapack_prebuild
    cd ..
fi
OPENBLAS_LIBCORE="$(grep 'LIBCORE=' $openblas_conf | cut -f2 -d=)"
OPENBLAS_ARCH="$(grep 'ARCH=' $openblas_conf | cut -f2 -d=)"
echo "OpenBLAS detected LIBCORE = $OPENBLAS_LIBCORE"
echo "OpenBLAS deteched ARCH    = $OPENBLAS_ARCH"
# output setup file
cat <<EOF > "${BUILDDIR}/openblas_arch"
export OPENBLAS_LIBCORE="${OPENBLAS_LIBCORE}"
export OPENBLAS_ARCH="${OPENBLAS_ARCH}"
EOF
