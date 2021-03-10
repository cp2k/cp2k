#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

mpich_ver="3.3.2"
mpich_sha256="4bfaf8837a54771d3e4922c84071ef80ffebddbb6971a006038d91ee7ef959b9"
mpich_pkg="mpich-${mpich_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ ${MPI_MODE} != "mpich" ] && exit 0
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
    if verify_checksums "${install_lock_file}"; then
      echo "mpich-${mpich_ver} is already installed, skipping it."
    else
      if [ -f ${mpich_pkg} ]; then
        echo "${mpich_pkg} is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${mpich_sha256} \
          https://www.cp2k.org/static/downloads/${mpich_pkg}
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d mpich-${mpich_ver} ] && rm -rf mpich-${mpich_ver}
      tar -xzf ${mpich_pkg}
      cd mpich-${mpich_ver}
      unset F90
      unset F90FLAGS

      # workaround for compilation with GCC-10, until properly fixed:
      #   https://github.com/pmodels/mpich/issues/4300
      ("${FC}" --version | grep -Eq 'GNU.+\s10\.') && compat_flag="-fallow-argument-mismatch" || compat_flag=""
      ./configure --prefix="${pkg_install_dir}" --libdir="${pkg_install_dir}/lib" MPICC="" FFLAGS="${FCFLAGS} ${compat_flag}" FCFLAGS="${FCFLAGS} ${compat_flag}" --without-x --enable-gl=no > configure.log 2>&1
      make -j $(get_nprocs) > make.log 2>&1
      make install > install.log 2>&1
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage1/$(basename ${SCRIPT_NAME})"
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
  __DONTUSE__) ;;

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
if [ "$with_mpich" != "__DONTUSE__" ]; then
  MPICH_LIBS="-lmpi -lmpicxx"
  if [ "$with_mpich" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_mpich"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_mpich" >> $SETUPFILE
    mpi_bin="$pkg_install_dir/bin/mpirun"
  else
    mpi_bin=mpirun
  fi
  # check MPICH version, versions less than 3.0 will get -D__MPI_VERSION=2 flag
  raw_version=$($mpi_bin --version |
    grep "Version:" | awk '{print $2}')
  major_version=$(echo $raw_version | cut -d '.' -f 1)
  minor_version=$(echo $raw_version | cut -d '.' -f 2)
  if [ $major_version -lt 3 ]; then
    mpi2_dflags="-D__MPI_VERSION=2"
  else
    mpi2_dflags=''
  fi
  cat << EOF >> "${BUILDDIR}/setup_mpich"
export MPI_MODE="${MPI_MODE}"
export MPICH_CFLAGS="${MPICH_CFLAGS}"
export MPICH_LDFLAGS="${MPICH_LDFLAGS}"
export MPICH_LIBS="${MPICH_LIBS}"
export MPI_CFLAGS="${MPICH_CFLAGS}"
export MPI_LDFLAGS="${MPICH_LDFLAGS}"
export MPI_LIBS="${MPICH_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__parallel ${mpi2_dflags}|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${MPICH_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${MPICH_LDFLAGS}|)"
export CP_LIBS="\${CP_LIBS} IF_MPI(${MPICH_LIBS}|)"
EOF
fi

load "${BUILDDIR}/setup_mpich"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "mpich"
