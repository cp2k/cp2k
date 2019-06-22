#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

libint_ver="2.5.0"

case "$LIBINT_LMAX" in
    4)
        libint_sha256="8b797a518fd5a0fa19c420c84794a7ec5225821ffc56d1666845a406d0b62982"
        ;;
    5)
        libint_sha256="5f2f1f1ee2cd04ec302415ea9b20e2ef01173e083d05b02cd633fd6f640e8e4c"
        ;;
    6)
        libint_sha256="9301591a97e6a91a26dcefdcdea7e36c946e310b5e9574a3d749b20700ad832b"
        ;;
    7)
        libint_sha256="d6ccce6de5964642093edb6eac139cc587aee826394d0c54d68e706ede05b2c2"
        ;;
    *)
       report_error "Unsupported value --libint-lmax=${LIBINT_LMAX}."
       exit 1
       ;;
esac

[ -f "${BUILDDIR}/setup_libint" ] && rm "${BUILDDIR}/setup_libint"

LIBINT_CFLAGS=''
LIBINT_LDFLAGS=''
LIBINT_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libint" in
    __INSTALL__)
        echo "==================== Installing LIBINT ===================="
        pkg_install_dir="${INSTALLDIR}/libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "libint-${libint_ver} is already installed, skipping it."
        else
            if [ -f libint_cp2k-${libint_ver}.tgz ] ; then
                echo "libint_cp2k-${libint_ver}.tgz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${libint_sha256} \
                             https://github.com/cp2k/libint-cp2k/releases/download/v${libint_ver}/libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}.tgz
            fi

            [ -d libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX} ] && rm -rf libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}
            tar -xzf libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}.tgz

            echo "Installing from scratch into ${pkg_install_dir}"
            cd libint-v${libint_ver}-cp2k-lmax-${LIBINT_LMAX}

            # reduce debug information to level 1 since
            # level 2 (default for -g flag) leads to very large binary size
            LIBINT_CXXFLAGS="$CXXFLAGS -g1"

            # hack for -with-cxx, needed for -fsanitize=thread that also
            # needs to be passed to the linker, but seemingly ldflags is
            # ignored by libint configure

            ./configure --prefix=${pkg_install_dir} \
                        --with-cxx="$CXX $LIBINT_CXXFLAGS" \
                        --with-cxx-optflags="$LIBINT_CXXFLAGS" \
                        --enable-fortran \
                        --libdir="${pkg_install_dir}/lib" \
                        > configure.log 2>&1

            make -j $NPROCS > make.log 2>&1
            make install > install.log 2>&1
            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi

        LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LIBINT from system paths ===================="
        check_lib -lint2 "libint"
        add_include_from_paths -p LIBINT_CFLAGS "libint" $INCLUDE_PATHS
        add_lib_from_paths LIBINT_LDFLAGS "libint2.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking LIBINT to user paths ===================="
        pkg_install_dir="$with_libint"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_libint" != "__DONTUSE__" ] ; then
    LIBINT_LIBS="-lint2"
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
export CP_DFLAGS="\${CP_DFLAGS} -D__LIBINT"
export CP_CFLAGS="\${CP_CFLAGS} ${LIBINT_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBINT_LDFLAGS}"
export CP_LIBS="${LIBINT_LIBS} \${CP_LIBS}"
EOF
fi

# update toolchain environment
load "${BUILDDIR}/setup_libint"
export -p > "${INSTALLDIR}/toolchain.env"

cd "${ROOTDIR}"
report_timing "libint"
