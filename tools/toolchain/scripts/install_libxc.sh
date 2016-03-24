#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_libxc=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_libxc" ] && rm "${BUILDDIR}/setup_libxc"

LIBXC_CFLAGS=''
LIBXC_LDFLAGS=''
LIBXC_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_libxc" in
    __INSTALL__)
        echo "==================== Installing LIBXC ===================="
        pkg_install_dir="${INSTALLDIR}/libxc-${libxc_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "libxc-${libxc_ver} is already installed, skipping it."
        else
            if [ -f libxc-${libxc_ver}.tar.gz ] ; then
                echo "libxc-${libxc_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/libxc-${libxc_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xzf libxc-${libxc_ver}.tar.gz
            cd libxc-${libxc_ver}
            # patch buggy configure macro (fails with gcc trunk)
            sed -i 's/ax_cv_f90_modext=$(ls | sed/ax_cv_f90_modext=)ls -1 | grep -iv smod | sed/g' \
                configure
            ./configure  --prefix="${pkg_install_dir}" --libdir="${pkg_install_dir}/lib" > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        LIBXC_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LIBXC from system paths ===================="
        check_lib -lxcf90 "libxc"
        check_lib -lxc "libxc"
        add_include_from_paths LIBXC_CFLAGS "xc.h" $INCLUDE_PATHS
        add_lib_from_paths LIBXC_LDFLAGS "libxc.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking LIBXC to user paths ===================="
        pkg_install_dir="$with_libxc"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        LIBXC_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBXC_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_libxc" != "__DONTUSE__" ] ; then
    LIBXC_LIBS="-lxcf90 -lxc"
    if [ "$with_libxc" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_libxc"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_libxc" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_libxc"
export LIBXC_CFLAGS="${LIBXC_CFLAGS}"
export LIBXC_LDFLAGS="${LIBXC_LDFLAGS}"
export LIBXC_LIBS="${LIBXC_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBXC"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBXC_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBXC_LDFLAGS}"
export CP_LIBS="${LIBXC_LIBS} \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
