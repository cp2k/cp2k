#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

mpich_ver="4.0.3"
mpich_sha256="17406ea90a6ed4ecd5be39c9ddcbfac9343e6ab4f77ac4e8c5ebe4a3e3b6c501"
mpich_pkg="mpich-${mpich_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "mpich" ] && exit 0
[ -f "${BUILDDIR}/setup_mpich" ] && rm "${BUILDDIR}/setup_mpich"

MPICH_CFLAGS=""
MPICH_LDFLAGS=""
MPICH_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_mpich}" in
  __INSTALL__)
    echo "==================== Installing MPICH ===================="
    pkg_install_dir="${INSTALLDIR}/mpich-${mpich_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "mpich-${mpich_ver} is already installed, skipping it."
    else
      if [ -f ${mpich_pkg} ]; then
        echo "${mpich_pkg} is found"
      else
        download_pkg_from_cp2k_org "${mpich_sha256}" "${mpich_pkg}"
      fi
      echo "Installing from scratch into ${pkg_install_dir} for MPICH device ${MPICH_DEVICE}"
      [ -d mpich-${mpich_ver} ] && rm -rf mpich-${mpich_ver}
      tar -xzf ${mpich_pkg}
      cd mpich-${mpich_ver}
      unset F90
      unset F90FLAGS

      # workaround for compilation with GCC >= 10, until properly fixed:
      #   https://github.com/pmodels/mpich/issues/4300
      if ("${FC}" --version | grep -q 'GNU'); then
        compat_flag=$(allowed_gfortran_flags "-fallow-argument-mismatch")
      fi
      ./configure \
        --prefix="${pkg_install_dir}" \
        --libdir="${pkg_install_dir}/lib" \
        MPICC="" \
        FFLAGS="${FCFLAGS} ${compat_flag}" \
        FCFLAGS="${FCFLAGS} ${compat_flag}" \
        --without-x --without-slurm \
        --enable-gl=no \
        --with-device=${MPICH_DEVICE} \
        > configure.log 2>&1 || tail -n ${LOG_LINES} configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      make install > install.log 2>&1 || tail -n ${LOG_LINES} install.log
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage1/$(basename ${SCRIPT_NAME})"
    fi
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_install ${pkg_install_dir}/bin/mpiexec "mpich" && MPIRUN="${pkg_install_dir}/bin/mpiexec" || exit 1
    check_install ${pkg_install_dir}/bin/mpicc "mpich" && MPICC="${pkg_install_dir}/bin/mpicc" || exit 1
    check_install ${pkg_install_dir}/bin/mpicxx "mpich" && MPICXX="${pkg_install_dir}/bin/mpicxx" || exit 1
    check_install ${pkg_install_dir}/bin/mpifort "mpich" && MPIFC="${pkg_install_dir}/bin/mpifort" || exit 1
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    MPICH_CFLAGS="-I'${pkg_install_dir}/include'"
    MPICH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding MPICH from system paths ===================="
    check_command mpiexec "mpich" && MPIRUN="$(command -v mpiexec)"
    check_command mpicc "mpich" && MPICC="$(command -v mpicc)" || exit 1
    if [ $(command -v mpic++ > /dev/null 2>&1) ]; then
      check_command mpic++ "mpich" && MPICXX="$(command -v mpic++)" || exit 1
    else
      check_command mpicxx "mpich" && MPICXX="$(command -v mpicxx)" || exit 1
    fi
    check_command mpifort "mpich" && MPIFC="$(command -v mpifort)" || exit 1
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    check_lib -lmpifort "mpich"
    check_lib -lmpicxx "mpich"
    check_lib -lmpi "mpich"
    add_include_from_paths MPICH_CFLAGS "mpi.h" ${INCLUDE_PATHS}
    add_lib_from_paths MPICH_LDFLAGS "libmpi.*" ${LIB_PATHS}
    ;;
  __DONTUSE__)
    # Nothing to do
    ;;
  *)
    echo "==================== Linking MPICH to user paths ===================="
    pkg_install_dir="${with_mpich}"
    check_dir "${pkg_install_dir}/bin"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    check_command ${pkg_install_dir}/bin/mpiexec "mpich" && MPIRUN="${pkg_install_dir}/bin/mpiexec" || exit 1
    check_command ${pkg_install_dir}/bin/mpicc "mpich" && MPICC="${pkg_install_dir}/bin/mpicc" || exit 1
    check_command ${pkg_install_dir}/bin/mpicxx "mpich" && MPICXX="${pkg_install_dir}/bin/mpicxx" || exit 1
    check_command ${pkg_install_dir}/bin/mpifort "mpich" && MPIFC="${pkg_install_dir}/bin/mpifort" || exit 1
    MPIFORT="${MPIFC}"
    MPIF77="${MPIFC}"
    MPICH_CFLAGS="-I'${pkg_install_dir}/include'"
    MPICH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "${with_mpich}" != "__DONTUSE__" ]; then
  if [ "${with_mpich}" != "__SYSTEM__" ]; then
    mpi_bin="${pkg_install_dir}/bin/mpiexec"
  else
    mpi_bin="mpiexec"
  fi
  MPICH_LIBS="-lmpifort -lmpicxx -lmpi"
  cat << EOF > "${BUILDDIR}/setup_mpich"
export MPI_MODE="${MPI_MODE}"
export MPIRUN="${MPIRUN}"
export MPICC="${MPICC}"
export MPICXX="${MPICXX}"
export MPIFC="${MPIFC}"
export MPIFORT="${MPIFORT}"
export MPIF77="${MPIF77}"
export MPICH_CFLAGS="${MPICH_CFLAGS}"
export MPICH_LDFLAGS="${MPICH_LDFLAGS}"
export MPICH_LIBS="${MPICH_LIBS}"
export MPI_CFLAGS="${MPICH_CFLAGS}"
export MPI_LDFLAGS="${MPICH_LDFLAGS}"
export MPI_LIBS="${MPICH_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${MPICH_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${MPICH_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_MPI(${MPICH_LIBS}|)"
EOF
  if [ "${with_mpich}" != "__SYSTEM__" ]; then
    # Using append_path instead of prepend_path for compatibility with Shifter.
    # See http://github.com/cp2k/cp2k/issues/2058
    cat << EOF >> "${BUILDDIR}/setup_mpich"
append_path PATH "${pkg_install_dir}/bin"
append_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
append_path LD_RUN_PATH "${pkg_install_dir}/lib"
append_path LIBRARY_PATH "${pkg_install_dir}/lib"
append_path CPATH "${pkg_install_dir}/include"
EOF
  fi
  cat "${BUILDDIR}/setup_mpich" >> ${SETUPFILE}
fi

# Update leak suppression file
cat << EOF >> ${INSTALLDIR}/lsan.supp
# MPICH 3.3.2 with GCC 10.3.0
leak:MPIR_Find_local_and_external
leak:MPIU_Find_local_and_external
EOF

load "${BUILDDIR}/setup_mpich"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "mpich"
