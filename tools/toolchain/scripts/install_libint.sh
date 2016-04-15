#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_libint=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_libint" ] && rm "${BUILDDIR}/setup_libint"

LIBINT_CFLAGS=''
LIBINT_LDFLAGS=''
LIBINT_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_libint" in
    __INSTALL__)
        echo "==================== Installing LIBINT ===================="
        pkg_install_dir="${INSTALLDIR}/libint-${libint_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "libint-${libint_ver} is already installed, skipping it."
        else
            if [ -f libint-${libint_ver}.tar.gz ] ; then
                echo "libint-${libint_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/libint-${libint_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d libint-${libint_ver} ] && rm -rf libint-${libint_ver}
            tar -xzf libint-${libint_ver}.tar.gz
            cd libint-${libint_ver}
            # hack for -with-cc, needed for -fsanitize=thread that also
            # needs to be passed to the linker, but seemingly ldflags is
            # ignored by libint configure
            ./configure --prefix=${pkg_install_dir} \
                        --libdir="${pkg_install_dir}/lib" \
                        --with-libint-max-am=5 \
                        --with-libderiv-max-am1=4 \
                        --with-cc="$CC $CFLAGS" \
                        --with-cc-optflags="$CFLAGS" \
                        --with-cxx-optflags="$CXXFLAGS" \
                        > configure.log 2>&1
            make -j $NPROCS >  make.log 2>&1
            make install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LIBINT from system paths ===================="
        check_lib -lderiv "libint"
        check_lib -lint "libint"
        add_include_from_paths -p LIBINT_CFLAGS "libint" $INCLUDE_PATHS
        add_lib_from_paths LIBINT_LDFLAGS "libint.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking LIBINT to user paths ===================="
        pkg_install_dir="$with_libint"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_libint" != "__DONTUSE__" ] ; then
    LIBINT_LIBS="-lderiv -lint"
    if [ "$with_libint" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_libint"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_libint" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_libint"
export LIBINT_CFLAGS="${LIBINT_CFLAGS}"
export LIBINT_LDFLAGS="${LIBINT_LDFLAGS}"
export LIBINT_LIBS="${LIBINT_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBINT -D__LIBINT_MAX_AM=6 -D__LIBDERIV_MAX_AM1=5"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBINT_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBINT_LDFLAGS}"
export CP_LIBS="${LIBINT_LIBS} \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
