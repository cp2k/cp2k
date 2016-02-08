#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

# helper to check if libsmm is available (uses https-redirect
# to find latest version)
libsmm_exists() {
    query_url=https://www.cp2k.org/static/downloads/libsmm/$1-latest.a
    reply_url="$(python -c "import urllib2; print(urllib2.urlopen('$query_url').geturl())" 2>&-)"
    if [ "$query_url" != "$reply_url" ]; then
        echo $reply_url | cut -d/ -f7
    fi
}

with_libsmm=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_libsmm" ] && rm "${BUILDDIR}/setup_libsmm"

LIBSMM_CFLAGS=''
LIBSMM_LDFLAGS=''
LIBSMM_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_libsmm" in
    __INSTALL__)
        echo "==================== Installing libsmm ===================="
        pkg_install_dir="${INSTALLDIR}/libsmm"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "libsmm is already installed, skipping it."
        else
            echo "Searching for an optimised libsmm binary from CP2K website"
            # Here we attempt to determine which precompiled libsmm binary
            # to download, and do that if such binary exists on CP2K web
            # repository.  The binary is determined by the arch and
            # libcore values obtained via OpenBLAS prebuild.
            libsmm="$(libsmm_exists libsmm_dnn_${OPENBLAS_LIBCORE})"
            if [ "x$libsmm" != "x" ] ; then
                echo "An optimized libsmm $libsmm is available"
            else
                echo "No optimised binary found ..."
                echo "Searching for a generic libsmm binary from CP2K website"
                libsmm="$(libsmm_exists libsmm_dnn_${OPENBLAS_ARCH})"
                if [ "x$libsmm" != "x" ] ; then
                    echo "A generic libsmm $libsmm is available."
                    echo "Consider building and contributing to CP2K an optimized"
                    echo "libsmm for your $OPENBLAS_ARCH $OPENBLAS_LIBCORE using"
                    echo "the toolkit in tools/build_libsmm provided in cp2k package"
                fi
            fi
            # we know what to get, proceed with install
            if [ "x$libsmm" != "x" ]; then
                if [ -f $libsmm ]; then
                    echo "$libsmm has already been downloaded."
                else
                    download_pkg ${DOWNLOADER_FLAGS} \
                                 https://www.cp2k.org/static/downloads/libsmm/$libsmm
                fi
                # install manually
                ! [ -d "${pkg_install_dir}/lib" ] && mkdir -p "${pkg_install_dir}/lib"
                cp $libsmm "${pkg_install_dir}/lib"
                ln -s "${pkg_install_dir}/lib/$libsmm" "${pkg_install_dir}/lib/libsmm_dnn.a"
            else
                echo "No libsmm is available"
                echo "Consider building an optimized libsmm on your system yourself"
                echo "using the toolkid in tools/build_libsmm provided in cp2k package"
                cat <<EOF > "${BUILDDIR}/setup_libsmm"
with_libsmm="__DONTUSE__"
EOF
                exit 0
            fi
            touch "${install_lock_file}"
        fi
        LIBSMM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding Libsmm from system paths ===================="
        check_lib -lsmm_dnn "libsmm"
        add_lib_from_paths LIBSMM_LDFLAGS "libsmm_dnn.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking Libsmm to user paths ===================="
        pkg_install_dir="$with_libsmm"
        check_dir "${pkg_install_dir}/lib"
        LIBSMM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_libsmm" != "__DONTUSE__" ] ; then
    LIBSMM_LIBS="-lsmm_dnn"
    if [ "$with_libsmm" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_libsmm"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
EOF
        cat "${BUILDDIR}/setup_libsmm" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_libsmm"
export LIBSMM_LDFLAGS="${LIBSMM_LDFLAGS}"
export LIBSMM_LIBS="${LIBSMM_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_VALGRIND(|-D__HAS_smm_dnn)"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBSMM_LDFLAGS}"
export CP_LIBS="IF_VALGRIND(|${LIBSMM_LIBS}) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
