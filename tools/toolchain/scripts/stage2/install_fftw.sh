#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

fftw_ver="3.3.8"
fftw_sha256="6113262f6e92c5bd474f2875fa1b01054c4ad5040f6b0da7c03c98821d9ae303"
fftw_pkg="fftw-${fftw_ver}.tar.gz"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_fftw" ] && rm "${BUILDDIR}/setup_fftw"

FFTW_CFLAGS=''
FFTW_LDFLAGS=''
FFTW_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_fftw" in
  __INSTALL__)
    require_env MPI_LIBS
    echo "==================== Installing FFTW ===================="
    pkg_install_dir="${INSTALLDIR}/fftw-${fftw_ver}"
    install_lock_file="$pkg_install_dir/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "fftw-${fftw_ver} is already installed, skipping it."
    else
      if [ -f ${fftw_pkg} ]; then
        echo "${fftw_pkg} is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${fftw_sha256} \
          "https://www.cp2k.org/static/downloads/${fftw_pkg}"
      fi
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d fftw-${fftw_ver} ] && rm -rf fftw-${fftw_ver}
      tar -xzf ${fftw_pkg}
      cd fftw-${fftw_ver}
      FFTW_FLAGS="--enable-openmp --enable-shared --enable-static"
      # fftw has mpi support but not compiled by default. so compile it if we build with mpi.
      # it will create a second library to link with if needed
      [ "$MPI_MODE" != "no" ] && FFTW_FLAGS="--enable-mpi ${FFTW_FLAGS}"
      grep '\bavx\b' /proc/cpuinfo 1> /dev/null && FFTW_FLAGS="${FFTW_FLAGS} --enable-avx"
      grep '\bavx2\b' /proc/cpuinfo 1> /dev/null && FFTW_FLAGS="${FFTW_FLAGS} --enable-avx2"
      grep '\bavx512f\b' /proc/cpuinfo 1> /dev/null && FFTW_FLAGS="${FFTW_FLAGS} --enable-avx512"
      ./configure --prefix=${pkg_install_dir} --libdir="${pkg_install_dir}/lib" ${FFTW_FLAGS} > configure.log 2>&1
      make -j $(get_nprocs) > make.log 2>&1
      make install > install.log 2>&1
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage2/$(basename ${SCRIPT_NAME})"
    fi
    FFTW_CFLAGS="-I'${pkg_install_dir}/include'"
    FFTW_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding FFTW from system paths ===================="
    check_lib -lfftw3 "FFTW"
    check_lib -lfftw3_omp "FFTW"
    [ "$MPI_MODE" != "no" ] && check_lib -lfftw3_mpi "FFTW"
    add_include_from_paths FFTW_CFLAGS "fftw3.h" $INCLUDE_PATHS
    add_lib_from_paths FFTW_LDFLAGS "libfftw3.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking FFTW to user paths ===================="
    pkg_install_dir="$with_fftw"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    FFTW_CFLAGS="-I'${pkg_install_dir}/include'"
    FFTW_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_fftw" != "__DONTUSE__" ]; then
  FFTW_LIBS="-lfftw3 -lfftw3_omp"
  if [ "$with_fftw" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_fftw"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
  fi
  # we may also want to cover FFT_SG
  cat << EOF >> "${BUILDDIR}/setup_fftw"
export FFTW_CFLAGS="${FFTW_CFLAGS}"
export FFTW_LDFLAGS="${FFTW_LDFLAGS}"
export FFTW_LIBS="${FFTW_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__FFTW3 IF_COVERAGE(IF_MPI(|-U__FFTW3)|)"
export CP_CFLAGS="\${CP_CFLAGS} ${FFTW_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${FFTW_LDFLAGS}"
export CP_LIBS="IF_MPI(-lfftw3_mpi|) ${FFTW_LIBS} \${CP_LIBS}"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
export FFTW_ROOT="$pkg_install_dir"
EOF

  cat "${BUILDDIR}/setup_fftw" >> $SETUPFILE
fi
cd "${ROOTDIR}"

# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat << EOF >> ${INSTALLDIR}/valgrind.supp
{
   <BuggyFFTW3>
   Memcheck:Addr32
   fun:cdot
   ...
   fun:invoke_solver
   fun:search0
}
EOF

load "${BUILDDIR}/setup_fftw"
write_toolchain_env "${INSTALLDIR}"

report_timing "fftw"
