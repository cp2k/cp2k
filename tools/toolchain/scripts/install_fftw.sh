#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

fftw_ver="3.3.8"
fftw_sha256="6113262f6e92c5bd474f2875fa1b01054c4ad5040f6b0da7c03c98821d9ae303"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_fftw" ] && rm "${BUILDDIR}/setup_fftw"

FFTW_CFLAGS=''
FFTW_LDFLAGS=''
FFTW_LIBS=''
FFTW_LIBS_OMP=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_fftw" in
    __INSTALL__)
        require_env MPI_LIBS
        echo "==================== Installing FFTW ===================="
        pkg_install_dir="${INSTALLDIR}/fftw-${fftw_ver}"
        install_lock_file="$pkg_install_dir/install_successful"

        if verify_checksums "${install_lock_file}" ; then
            echo "fftw-${fftw_ver} is already installed, skipping it."
        else
            if [ -f fftw-${fftw_ver}.tar.gz ] ; then
                echo "fftw-${fftw_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${fftw_sha256} \
                             "https://www.cp2k.org/static/downloads/fftw-${fftw_ver}.tar.gz"
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d fftw-${fftw_ver} ] && rm -rf fftw-${fftw_ver}
            tar -xzf fftw-${fftw_ver}.tar.gz
            cd fftw-${fftw_ver}
            if [ "$MPI_MODE" != "no" ] ; then
                # fftw has mpi support but not compiled by default. so compile it if we build with mpi.
                # it will create a second library to link with if needed
                ./configure  --prefix=${pkg_install_dir} --libdir="${pkg_install_dir}/lib" --enable-openmp --enable-mpi --enable-shared > configure.log 2>&1
            else
                ./configure  --prefix=${pkg_install_dir} --libdir="${pkg_install_dir}/lib" --enable-openmp --enable-shared > configure.log 2>&1
            fi
            make -j $NPROCS > make.log 2>&1
            make install > install.log 2>&1
            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi
        FFTW_CFLAGS="-I'${pkg_install_dir}/include'"
        FFTW_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding FFTW from system paths ===================="
        check_lib -lfftw3 "FFTW"
        [ $ENABLE_OMP = "__TRUE__" ] && check_lib -lfftw3_omp "FFTW"
        [ "$MPI_MODE" != "no" ] && check_lib -lfftw3_mpi "FFTW"
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
export CP_LIBS="IF_MPI(-lfftw3_mpi|) ${FFTW_LIBS} IF_OMP(${FFTW_LIBS_OMP}|) \${CP_LIBS}"
prepend_path PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$pkg_install_dir/lib/pkgconfig"
export FFTWROOT="$pkg_install_dir"
EOF
fi
cd "${ROOTDIR}"


# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat <<EOF >> ${INSTALLDIR}/valgrind.supp
{
   <BuggyFFTW3>
   Memcheck:Addr32
   fun:cdot
   ...
   fun:invoke_solver
   fun:search0
}
EOF

# update toolchain environment
load "${BUILDDIR}/setup_fftw"
export -p > "${INSTALLDIR}/toolchain.env"

report_timing "fftw"
