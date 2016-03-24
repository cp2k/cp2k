#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_fftw=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_fftw" ] && rm "${BUILDDIR}/setup_fftw"

FFTW_CFLAGS=''
FFTW_LDFLAGS=''
FFTW_LIBS=''
FFTW_LIBS_OMP=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_fftw" in
    __INSTALL__)
        echo "==================== Installing FFTW ===================="
        pkg_install_dir="${INSTALLDIR}/fftw-${fftw_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "fftw-${fftw_ver} is already installed, skipping it."
        else
            if [ -f fftw-${fftw_ver}.tar.gz ] ; then
                echo "fftw-${fftw_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/fftw-${fftw_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xzf fftw-${fftw_ver}.tar.gz
            cd fftw-${fftw_ver}
            ./configure  --prefix=${pkg_install_dir} --libdir="${pkg_install_dir}/lib" --enable-openmp > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        FFTW_CFLAGS="-I'${pkg_install_dir}/include'"
        FFTW_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding FFTW from system paths ===================="
        check_lib -lfftw3 "FFTW"
        [ $ENABLE_OMP = "__TRUE__" ] && check_lib -lfftw3_omp "FFTW"
        add_include_from_paths FFTW_CFLAGS "fftw3.h" $INCLUDE_PATHS
        add_lib_from_paths FFTW_LDFLAGS "libfftw3.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking FFTW to user paths ===================="
        pkg_install_dir="$with_fftw"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        FFTW_CFLAGS="-I'${pkg_install_dir}/include'"
        FFTW_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_fftw" != "__DONTUSE__" ] ; then
    FFTW_LIBS="-lfftw3"
    FFTW_LIBS_OMP="-lfftw3_omp"
    if [ "$with_fftw" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_fftw"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_fftw" >> $SETUPFILE
    fi
    # we may also want to cover FFT_SG
    cat <<EOF >> "${BUILDDIR}/setup_fftw"
export FFTW_CFLAGS="${FFTW_CFLAGS}"
export FFTW_LDFLAGS="${FFTW_LDFLAGS}"
export FFTW_LIBS="${FFTW_LIBS}"
export FFTW_LIBS_OMP="${FFTW_LIBS_OMP}"
export CP_DFLAGS="\${CP_DFLAGS} -D__FFTW3 IF_COVERAGE(IF_MPI(|-U__FFTW3)|)"
export CP_CFLAGS="\${CP_CFLAGS} ${FFTW_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${FFTW_LDFLAGS}"
export CP_LIBS="${FFTW_LIBS} IF_OMP(${FFTW_LIBS_OMP}|) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
