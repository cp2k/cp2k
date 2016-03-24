#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_openmpi=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_openmpi" ] && rm "${BUILDDIR}/setup_openmpi"

OPENMPI_CFLAGS=''
OPENMPI_LDFLAGS=''
OPENMPI_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_openmpi" in
    __INSTALL__)
        echo "==================== Installing OpenMPI ===================="
        pkg_install_dir="${INSTALLDIR}/openmpi-${openmpi_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "openmpi-${openmpi_ver} is already installed, skipping it."
        else
            if [ -f openmpi-${openmpi_ver}.tar.gz ] ; then
                echo "openmpi-${openmpi_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/openmpi-${openmpi_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xzf openmpi-${openmpi_ver}.tar.gz
            cd openmpi-${openmpi_ver}
            # can have issue with older glibc libraries, in which case
            # we need to add the -fgnu89-inline to CFLAGS. We can check
            # the version of glibc using ldd --version, as ldd is part of
            # glibc package
            glibc_version=$(ldd --version | awk '(NR == 1){print $4}')
            glibc_major_ver=$(echo $glibc_version | cut -d . -f 1)
            glibc_minor_ver=$(echo $glibc_version | cut -d . -f 2)
            if [ $glibc_major_ver -lt 2 ] || \
               [ $glibc_major_ver -eq 2 -a $glibc_minor_ver -lt 12 ] ; then
                CFLAGS="${CFLAGS} -fgnu89-inline"
            fi
            ./configure --prefix=${pkg_install_dir} --libdir="${pkg_install_dir}/lib" CFLAGS="${CFLAGS}" > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        OPENMPI_CFLAGS="-I'${pkg_install_dir}/include'"
        OPENMPI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding OpenMPI from system paths ===================="
        check_command mpirun "openmpi"
        check_command mpicc "openmpi"
        check_command mpif90 "openmpi"
        check_command mpic++ "openmpi"
        check_lib -lmpi "openmpi"
        check_lib -lmpi_cxx "openmpi"
        add_include_from_paths OPENMPI_CFLAGS "mpi.h" $INCLUDE_PATHS
        add_lib_from_paths OPENMPI_LDFLAGS "libmpi.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking OpenMPI to user paths ===================="
        pkg_install_dir="$with_openmpi"
        check_dir "${pkg_install_dir}/bin"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        OPENMPI_CFLAGS="-I'${pkg_install_dir}/include'"
        OPENMPI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_openmpi" != "__DONTUSE__" ] ; then
    OPENMPI_LIBS="-lmpi -lmpi_cxx"
    if [ "$with_openmpi" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_openmpi"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_openmpi" >> $SETUPFILE
        mpi_bin="$pkg_install_dir/bin/mpirun"
    else
        mpi_bin=mpirun
    fi
    # check openmpi version, versions less than 1.7 will get -D__MPI_VERSION=2 flag
    raw_version=$($mpi_bin --version 2>&1 | \
                      grep "(Open MPI)" | awk '{print $4}')
    major_version=$(echo $raw_version | cut -d '.' -f 1)
    minor_version=$(echo $raw_version | cut -d '.' -f 2)
    if [ $major_version -lt 1 ] || \
       [ $major_version -eq 1 -a $minor_version -lt 7 ] ; then
        mpi2_dflags="-D__MPI_VERSION=2"
    else
        mpi2_dflags=''
    fi
    cat <<EOF >> "${BUILDDIR}/setup_openmpi"
export OPENMPI_CFLAGS="${OPENMPI_CFLAGS}"
export OPENMPI_LDFLAGS="${OPENMPI_LDFLAGS}"
export OPENMPI_LIBS="${OPENMPI_LIBS}"
export MPI_CFLAGS="${OPENMPI_CFLAGS}"
export MPI_LDFLAGS="${OPENMPI_LDFLAGS}"
export MPI_LIBS="${OPENMPI_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel ${mpi2_dflags}|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${OPENMPI_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${OPENMPI_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_MPI(${OPENMPI_LIBS}|)"
EOF
fi
cd "${ROOTDIR}"
