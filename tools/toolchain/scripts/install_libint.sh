#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

libint_ver="2.5.0"
libint_sha256="e57bb4546a6702fdaa570ad6607712f31903ed4618f051150979a31a038ce960"
boost_ver="1_70_0"
boost_sha256="882b48708d211a5f48e60b0124cf5863c1534cd544ecd0664bb534a4b5d506e9"
gmp_ver="6.1.2"
gmp_sha256="5275bb04f4863a13516b2f39392ac5e272f5e1bb8057b18aec1c9b79d73d8fb2"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

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
        if verify_checksums "${install_lock_file}" ; then
            echo "libint-${libint_ver} is already installed, skipping it."
        else
            if [ -f v${libint_ver}.tar.gz ] ; then
                echo "v${libint_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${libint_sha256} \
                             https://github.com/evaleev/libint/archive/v${libint_ver}.tar.gz
            fi
            if [ -f boost_${boost_ver}.tar.gz ] ; then
                echo "boost_${boost_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${boost_sha256} \
                             https://dl.bintray.com/boostorg/release/1.70.0/source/boost_${boost_ver}.tar.gz
            fi

            if [ -f gmp-${gmp_ver}.tar.bz2 ] ; then
                echo "gmp-${gmp_ver}.tar.bz2 is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${gmp_sha256} \
                             https://ftp.gnu.org/gnu/gmp/gmp-${gmp_ver}.tar.bz2
            fi

            [ -d libint-${libint_ver} ] && rm -rf libint-${libint_ver}
            tar -xzf v${libint_ver}.tar.gz
            [ -d boost_${boost_ver} ] && rm -rf boost_${boost_ver}
            tar -xzf boost_${boost_ver}.tar.gz
            [ -d gmp-${gmp_ver} ] && rm -rf gmp-${gmp_ver}
            tar -xjf gmp-${gmp_ver}.tar.bz2

            echo "Building dependency gmp-${gmp_ver}"

            cd gmp-${gmp_ver}
            ./configure --enable-cxx > configure.log 2>&1
            make > make.log 2>&1
            cd ..

            echo "Installing from scratch into ${pkg_install_dir}"
            cd libint-${libint_ver}
            ./autogen.sh
            mkdir build; cd build
            # hack for -with-cxx, needed for -fsanitize=thread that also
            # needs to be passed to the linker, but seemingly ldflags is
            # ignored by libint configure
            ../configure --enable-eri=1 --enable-eri3=1 \
                         --with-max-am=4 --with-eri-max-am=5,4 --with-eri3-max-am=5,4 \
                         --with-opt-am=3 --enable-generic-code --disable-unrolling \
                         --with-boost-libdir="${BUILDDIR}/boost_${boost_ver}/boost" \
                         --with-cxx="$CXX $CXXFLAGS" \
                         --with-cxx-optflags="$CXXFLAGS" \
                         --with-cxxgen-optflags="$CXXFLAGS" \
                         --with-incdirs="-I${BUILDDIR}/gmp-${gmp_ver}" \
                         --with-libdirs="-L${BUILDDIR}/gmp-${gmp_ver}/.libs"
                         #--with-libs="-lgmpxx -lgmp"
                         #> configure.log 2>&1

            make -j $NPROCS export > make.log 2>&1

            tar -xzf libint-${libint_ver}.tgz
            cd libint-${libint_ver}
            ./configure --prefix=${pkg_install_dir} \
                        --with-cxx="$CXX $CXXFLAGS" \
                        --with-cxx-optflags="$CXXFLAGS" \
                        --enable-fortran > configure.log 2>&1

            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS fortran >> make.log 2>&1
            make install > install.log 2>&1
            cd ../../..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi

        LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib64'"
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
        check_dir "${pkg_install_dir}/lib64"
        check_dir "${pkg_install_dir}/include"
        LIBINT_CFLAGS="-I'${pkg_install_dir}/include'"
        LIBINT_LDFLAGS="-L'${pkg_install_dir}/lib64'"
        ;;
esac
if [ "$with_libint" != "__DONTUSE__" ] ; then
    LIBINT_LIBS="-lint2"
    if [ "$with_libint" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_libint"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib64"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib64"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib64"
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
