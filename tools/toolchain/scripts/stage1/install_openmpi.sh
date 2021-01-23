#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

openmpi_ver="4.0.5"
openmpi_sha256="572e777441fd47d7f06f1b8a166e7f44b8ea01b8b2e79d1e299d509725d1bd05"
openmpi_pkg="openmpi-${openmpi_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "openmpi" ] && exit 0
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
        if verify_checksums "${install_lock_file}" ; then
            echo "openmpi-${openmpi_ver} is already installed, skipping it."
        else
            if [ -f ${openmpi_pkg} ] ; then
                echo "${openmpi_pkg} is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${openmpi_sha256} \
                             "https://www.cp2k.org/static/downloads/${openmpi_pkg}"
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d openmpi-${openmpi_ver} ] && rm -rf openmpi-${openmpi_ver}
            tar -xzf ${openmpi_pkg}
            cd openmpi-${openmpi_ver}
            # can have issue with older glibc libraries, in which case
            # we need to add the -fgnu89-inline to CFLAGS. We can check
            # the version of glibc using ldd --version, as ldd is part of
            # glibc package
            glibc_version=$(ldd --version | awk '/ldd/{print $NF}')
            glibc_major_ver=${glibc_version%%.*}
            glibc_minor_ver=${glibc_version##*.}
            if [ $glibc_major_ver -lt 2 ] || \
               [ $glibc_major_ver -eq 2 -a $glibc_minor_ver -lt 12 ] ; then
                CFLAGS="${CFLAGS} -fgnu89-inline"
            fi
            ./configure --prefix=${pkg_install_dir} --libdir="${pkg_install_dir}/lib" --enable-mpi1-compatibility CFLAGS="${CFLAGS}" > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage1/$(basename ${SCRIPT_NAME})"
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
        # Fortran code in CP2K is built via the mpifort wrapper, but we may need additional
        # libraries and linker flags for C/C++-based MPI codepaths, pull them in at this point.
        OPENMPI_CFLAGS="$(mpicxx --showme:compile)"
        OPENMPI_LDFLAGS="$(mpicxx --showme:link)"
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
        mpicxx_bin="$pkg_install_dir/bin/mpicxx"
    else
        mpi_bin=mpirun
        mpicxx_bin=mpicxx
    fi
    # check openmpi version as reported by mpirun
    raw_version=$($mpi_bin --version 2>&1 | \
                      grep "(Open MPI)" | awk '{print $4}')
    major_version=$(echo $raw_version | cut -d '.' -f 1)
    minor_version=$(echo $raw_version | cut -d '.' -f 2)
    OPENMPI_LIBS=""
    # grab additional runtime libs (for C/C++) from the mpicxx wrapper,
    # and remove them from the LDFLAGS if present
    for lib in $("${mpicxx_bin}" --showme:libs) ; do
        OPENMPI_LIBS+=" -l${lib}"
        OPENMPI_LDFLAGS="${OPENMPI_LDFLAGS//-l${lib}}"
    done
    # old versions didn't support MPI 3, so adjust __MPI_VERSION accordingly (needed e.g. for pexsi)
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


# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat <<EOF >> ${INSTALLDIR}/valgrind.supp
{
   <BuggyOpenMPI_1>
   Memcheck:Leak
   ...
   fun:*alloc
   ...
   fun:ompi_mpi_init
}
{
   <BuggyOpenMPI_2>
   Memcheck:Leak
   ...
   fun:*alloc
   ...
   fun:ompi_mpi_finalize
}
{
   <BuggyOpenMPI_3>
   Memcheck:Leak
   ...
   fun:malloc
   fun:opal_free_list_grow_st
   ...
   fun:mpi_alloc_mem
}
{
   <BuggyOpenMPI_4>
   Memcheck:Leak
   ...
   fun:malloc
   ...
   fun:progress_engine
   ...
   fun:clone
}
{
   <BuggyOpenMPI_5>
   Memcheck:Leak
   ...
   fun:malloc
   ...
   fun:query_2_0_0
   ...
   fun:ompi_comm_activate
}
EOF
cat <<EOF >> ${INSTALLDIR}/lsan.supp
# leaks related to OpenMPI
leak:query_2_0_0
leak:ompi_init_f
leak:ompi_finalize_f
leak:ompi_file_open_f
leak:progress_engine
leak:__GI___strdup
EOF

load "${BUILDDIR}/setup_openmpi"
write_toolchain_env "${INSTALLDIR}"

report_timing "openmpi"
