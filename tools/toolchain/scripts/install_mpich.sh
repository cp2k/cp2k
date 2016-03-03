#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_mpich=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_mpich" ] && rm "${BUILDDIR}/setup_mpich"

MPICH_CFLAGS=''
MPICH_LDFLAGS=''
MPICH_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_mpich" in
    __INSTALL__)
        echo "==================== Installing MPICH ===================="
        pkg_install_dir="${INSTALLDIR}/mpich-${mpich_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "mpich-${mpich_ver} is already installed, skipping it."
        else
            if [ -f mpich-${mpich_ver}.tar.gz ] ; then
                echo "mpich-${mpich_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/mpich-${mpich_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xzf mpich-${mpich_ver}.tar.gz
            cd mpich-${mpich_ver}
            (
                unset F90
                unset F90FLAGS
                ./configure --prefix="${pkg_install_dir}" > configure.log 2>&1
                make -j $NPROCS > make.log 2>&1
                make -j $NPROCS install > install.log 2>&1
            )
            cd ..
            touch "${install_lock_file}"
        fi
        MPICH_CFLAGS="-I'${pkg_install_dir}/include'"
        MPICH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding MPICH from system paths ===================="
        check_command mpirun "mpich"
        check_command mpicc "mpich"
        check_command mpif90 "mpich"
        check_command mpic++ "mpich"
        check_lib -lmpi "mpich"
        check_lib -lmpicxx "mpich"
        add_include_from_paths MPICH_CFLAGS "mpi.h" $INCLUDE_PATHS
        add_lib_from_paths MPICH_LDFLAGS "libmpi.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking MPICH to user paths ===================="
        pkg_install_dir="$with_mpich"
        check_dir "${pkg_install_dir}/bin"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        MPICH_CFLAGS="-I'${pkg_install_dir}/include'"
        MPICH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_mpich" != "__DONTUSE__" ] ; then
    MPICH_LIBS="-lmpi -lmpicxx"
    if [ "$with_mpich" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_mpich"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_mpich" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_mpich"
export MPI_MODE="__MPICH__"
export MPICH_CFLAGS="${MPICH_CFLAGS}"
export MPICH_LDFLAGS="${MPICH_LDFLAGS}"
export MPICH_LIBS="${MPICH_LIBS}"
export MPI_CFLAGS="${MPICH_CFLAGS}"
export MPI_LDFLAGS="${MPICH_LDFLAGS}"
export MPI_LIBS="${MPICH_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel -D__MPI_VERSION=3|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${MPICH_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${MPICH_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_MPI(${MPICH_LIBS}|)"
EOF
fi
cd "${ROOTDIR}"
