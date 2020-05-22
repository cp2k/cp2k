#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

gsl_ver="2.5"
gls_sha256="0460ad7c2542caaddc6729762952d345374784100223995eb14d614861f2258d"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_gsl" ] && rm "${BUILDDIR}/setup_gsl"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_gsl" in
    __INSTALL__)
        echo "==================== Installing gsl ===================="
        pkg_install_dir="${INSTALLDIR}/gsl-${gsl_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "gsl-${gsl_ver} is already installed, skipping it."
        else
            if [ -f gsl-${gsl_ver}.tar.gz ] ; then
                echo "gsl-${gsl_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${gls_sha256} \
                             "https://www.cp2k.org/static/downloads/gsl-${gsl_ver}.tar.gz"
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d gsl-${gsl_ver} ] && rm -rf gsl-${gsl_ver}
            tar -xzf gsl-${gsl_ver}.tar.gz
            cd gsl-${gsl_ver}
            ./configure --prefix="${pkg_install_dir}" \
                        --libdir="${pkg_install_dir}/lib" \
                        --enable-shared \
                        --enable-static > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi

        GSL_CFLAGS="-I'${pkg_install_dir}/include'"
        GSL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding gsl from system paths ===================="
        check_command pkg-config --modversion gsl
        add_include_from_paths GSL_CFLAGS "gsl.h" $INCLUDE_PATHS
        add_lib_from_paths GSL_LDFLAGS "libgsl.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking gsl to user paths ===================="
        pkg_install_dir="$with_gsl"
        check_dir "$pkg_install_dir/lib"
        check_dir "$pkg_install_dir/include"
        GSL_CFLAGS="-I'${pkg_install_dir}/include'"
        GSL_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_gsl" != "__DONTUSE__" ] ; then
    GSL_LIBS="-lgsl -lgslcblas"
    if [ "$with_gsl" != "__SYSTEM__" ] ; then
        cat << EOF > "${BUILDDIR}/setup_gsl"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
export GSL_INCLUDE_DIR="$pkg_install_dir/include"
export GSL_LIBRARY="-lgsl"
export GSL_CBLAS_LIBRARY="-lgslcblas"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$pkg_install_dir/lib64/pkgconfig:$pkg_install_dir/lib/pkgconfig"
EOF
    fi
    cat << EOF >> "${BUILDDIR}/setup_gsl"
export GSL_CFLAGS="${GSL_CFLAGS}"
export GSL_LDFLAGS="${GSL_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__GSL|)"
export CP_CFLAGS="\${CP_CFLAGS} ${GSL_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${GSL_LDFLAGS}"
export GSL_LIBRARY="-lgsl"
export GSLROOT="$pkg_install_dir"
export GSL_CBLAS_LIBRARY="-lgslcblas"
export GSL_INCLUDE_DIR="$pkg_install_dir/include"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$pkg_install_dir/lib64/pkgconfig:$pkg_install_dir/lib/pkgconfig"

##########################################################
#
# I only include the library when SIRIUS is activated
# which depends explicitly on MPI
#
##########################################################


export CP_LIBS="IF_MPI(${GSL_LIBS}|) \${CP_LIBS}"
EOF
    cat "${BUILDDIR}/setup_gsl" >> $SETUPFILE
fi

# update toolchain environment
load "${BUILDDIR}/setup_gsl"
export -p > "${INSTALLDIR}/toolchain.env"

cd "${ROOTDIR}"
report_timing "gsl"
