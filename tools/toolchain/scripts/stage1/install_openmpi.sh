#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

openmpi_ver="5.0.10"
openmpi_sha256="0acecc4fc218e5debdbcb8a41d182c6b0f1d29393015ed763b2a91d5d7374cc6"
openmpi_pkg="openmpi-${openmpi_ver}.tar.bz2"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "openmpi" ] && exit 0
[ -f "${BUILDDIR}/setup_openmpi" ] && rm "${BUILDDIR}/setup_openmpi"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_openmpi}" in
  __INSTALL__)
    echo "==================== Installing OpenMPI ===================="
    pkg_install_dir="${INSTALLDIR}/openmpi-${openmpi_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "openmpi-${openmpi_ver} is already installed, skipping it."
    else
      retrieve_package "${openmpi_sha256}" "${openmpi_pkg}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d openmpi-${openmpi_ver} ] && rm -rf openmpi-${openmpi_ver}
      tar -xjf ${openmpi_pkg}
      cd openmpi-${openmpi_ver}
      if [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
        # can have issue with older glibc libraries, in which case
        # we need to add the -fgnu89-inline to CFLAGS. We can check
        # the version of glibc using ldd --version, as ldd is part of
        # glibc package
        glibc_version=$(ldd --version | awk '/ldd/{print $NF}')
        glibc_major_ver=${glibc_version%%.*}
        glibc_minor_ver=${glibc_version##*.}
        if [ $glibc_major_ver -lt 2 ] ||
          [ $glibc_major_ver -eq 2 -a $glibc_minor_ver -lt 12 ]; then
          CFLAGS="${CFLAGS} -fgnu89-inline"
        fi
      fi

      ./configure CFLAGS="${CFLAGS}" \
        --prefix=${pkg_install_dir} \
        --libdir="${pkg_install_dir}/lib" \
        --enable-static \
        --with-libevent=internal \
        > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make -j $(get_nprocs) install > install.log 2>&1 || tail_excerpt install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage1/$(basename ${SCRIPT_NAME})"
    fi
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_install ${pkg_install_dir}/bin/mpiexec "openmpi" && MPIEXEC="${pkg_install_dir}/bin/mpiexec" || exit 1
    check_install ${pkg_install_dir}/bin/mpicc "openmpi" && MPICC="${pkg_install_dir}/bin/mpicc" || exit 1
    check_install ${pkg_install_dir}/bin/mpicxx "openmpi" && MPICXX="${pkg_install_dir}/bin/mpicxx" || exit 1
    check_install ${pkg_install_dir}/bin/mpifort "openmpi" && MPIFC="${pkg_install_dir}/bin/mpifort" || exit 1
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    ;;
  __SYSTEM__)
    echo "==================== Finding OpenMPI from system paths ===================="
    check_command mpiexec "openmpi" && MPIEXEC="$(command -v mpiexec)" || exit 1
    check_command mpicc "openmpi" && MPICC="$(command -v mpicc)" || exit 1
    check_command mpic++ "openmpi" && MPICXX="$(command -v mpic++)" || exit 1
    check_command mpifort "openmpi" && MPIFC="$(command -v mpifort)" || exit 1
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking OpenMPI to user paths ===================="
    pkg_install_dir="${with_openmpi}"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_command ${pkg_install_dir}/bin/mpiexec "openmpi" && MPIEXEC="${pkg_install_dir}/bin/mpiexec" || exit 1
    check_command ${pkg_install_dir}/bin/mpicc "openmpi" && MPICC="${pkg_install_dir}/bin/mpicc" || exit 1
    check_command ${pkg_install_dir}/bin/mpic++ "openmpi" && MPICXX="${pkg_install_dir}/bin/mpic++" || exit 1
    check_command ${pkg_install_dir}/bin/mpifort "openmpi" && MPIFC="${pkg_install_dir}/bin/mpifort" || exit 1
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    ;;
esac
if [ "${with_openmpi}" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_openmpi"
export MPI_MODE="${MPI_MODE}"
export MPIEXEC="${MPIEXEC}"
export MPICC="${MPICC}"
export MPICXX="${MPICXX}"
export MPIFC="${MPIFC}"
export MPIFORT="${MPIFORT}"
export MPIF77="${MPIF77}"
EOF
  if [ "${with_openmpi}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_openmpi"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_openmpi" "${SETUPFILE}"
fi

# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat << EOF >> ${INSTALLDIR}/valgrind.supp
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
cat << EOF >> ${INSTALLDIR}/lsan.supp
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

cd "${ROOTDIR}"
report_timing "openmpi"
