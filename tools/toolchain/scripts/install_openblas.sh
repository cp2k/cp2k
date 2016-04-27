#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_openblas=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_openblas" ] && rm "${BUILDDIR}/setup_openblas"

OPENBLAS_CFLAGS=''
OPENBLAS_LDFLAGS=''
OPENBLAS_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_openblas" in
    __INSTALL__)
        echo "==================== Installing OpenBLAS ===================="
        pkg_install_dir="${INSTALLDIR}/openblas-${openblas_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "openblas-${openblas_ver} is already installed, skipping it."
        else
            if [ -f OpenBLAS-${openblas_ver}.tar.gz ] ; then
                echo "OpenBLAS-${openblas_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/OpenBLAS-${openblas_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d OpenBLAS-${openblas_ver} ] && rm -rf OpenBLAS-${openblas_ver}
            tar -zxf OpenBLAS-${openblas_ver}.tar.gz
            cd OpenBLAS-${openblas_ver}
            # Originally we try to install both the serial and the omp
            # threaded version. Unfortunately, neither is thread-safe
            # (i.e. the CP2K ssmp and psmp version need to link to
            # something else, the omp version is unused)
            #
            # First attempt to make openblas using auto detected
            # TARGET, if this fails, then make with forced
            # TARGET=NEHALEM
            #
            # basename : see https://github.com/xianyi/OpenBLAS/issues/857
            #
            ( make -j $NPROCS \
                   USE_THREAD=0 \
                   CC=$(basename $CC) \
                   FC=$(basename $FC) \
                   PREFIX="${pkg_install_dir}" \
                   > make.serial.log 2>&1 \
            ) || ( \
                make -j $NPROCS clean; \
                make -j $NPROCS \
                     TARGET=NEHALEM \
                     USE_THREAD=0 \
                     CC=$(basename $CC) \
                     FC=$(basename $FC) \
                     PREFIX="${pkg_install_dir}" \
                     > make.serial.log 2>&1 \
            )
            make -j $NPROCS \
                 USE_THREAD=0 \
                 CC=$(basename $CC) \
                 FC=$(basename $FC) \
                 PREFIX="${pkg_install_dir}" \
                 install > install.serial.log 2>&1
            # make clean > clean.log 2>&1
            # make -j $nprocs \
            #      USE_THREAD=1 \
            #      USE_OPENMP=1 \
            #      LIBNAMESUFFIX=omp \
            #      CC=$(basename $CC) \
            #      FC=$(basename $FC) \
            #      PREFIX="${pkg_install_dir}" \
            #      > make.omp.log 2>&1
            # make -j $nprocs \
            #      USE_THREAD=1 \
            #      USE_OPENMP=1 \
            #      LIBNAMESUFFIX=omp \
            #      CC=$(basename $CC) \
            #      FC=$(basename $FC) \
            #      PREFIX="${pkg_install_dir}" \
            #      install > install.omp.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        OPENBLAS_CFLAGS="-I'${pkg_install_dir}/include'"
        OPENBLAS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding LAPACK from system paths ===================="
        check_lib -lopenblas "OpenBLAS"
        add_include_from_paths OPENBLAS_CFLAGS "openblas_config.h" $INCLUDE_PATHS
        add_lib_from_paths OPENBLAS_LDFLAGS "libopenblas.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking LAPACK to user paths ===================="
        pkg_install_dir="$with_openblas"
        check_dir "${pkg_install_dir}/include"
        check_dir "${pkg_install_dir}/lib"
        OPENBLAS_CFLAGS="-I'${pkg_install_dir}/include'"
        OPENBLAS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_openblas" != "__DONTUSE__" ] ; then
    OPENBLAS_LIBS="-lopenblas"
    if [ "$with_openblas" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_openblas"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_openblas" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_openblas"
export OPENBLAS_CFLAGS="${OPENBLAS_CFLAGS}"
export OPENBLAS_LDFLAGS="${OPENBLAS_LDFLAGS}"
export OPENBLAS_LIBS="${OPENBLAS_LIBS}"
export FAST_MATH_CFLAGS="\${FAST_MATH_CFLAGS} ${OPENBLAS_CFLAGS}"
export FAST_MATH_LDFLAGS="\${FAST_MATH_LDFLAGS} ${OPENBLAS_LDFLAGS}"
export FAST_MATH_LIBS="\${FAST_MATH_LIBS} ${OPENBLAS_LIBS}"
EOF
fi
cd "${ROOTDIR}"
