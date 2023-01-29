#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

DBCSR_ver="2.5.0"
DBCSR_sha256="e5c545ec16688027537f7865976b905c0783d038ec289e65635e63e961330601"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_dbcsr" ] && rm "${BUILDDIR}/setup_dbcsr"

DBCSR_CFLAGS=''
DBCSR_LDFLAGS=''
DBCSR_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

####case "$with_dbcsr" in
####    __INSTALL__)
echo "==================== Installing DBCSR ===================="
#
# to be restored to the right value when this script is included in the toolchain
pkg_install_dir="${INSTALLDIR}/DBCSR"
install_lock_file="$pkg_install_dir/install_successful"
if verify_checksums "${install_lock_file}"; then
  echo "DBCSR-${DBCSR_ver} is already installed, skipping it."
else
  if [ -f dbcsr-${DBCSR_ver}.tar.gz ]; then
    echo "dbcsr-${DBCSR_ver}.tar.gz is found"
  else
    wget -q "https://github.com/cp2k/dbcsr/archive/refs/tags/v${DBCSR_ver}.tar.gz" -O "dbcsr-${DBCSR_ver}.tar.gz"
  fi
  echo "Installing from scratch into ${pkg_install_dir}"
  [ -d dbcsr-${DBCSR_ver} ] && rm -rf dbcsr-${DBCSR_ver}
  tar xzf dbcsr-${DBCSR_ver}.tar.gz
  cd dbcsr-${DBCSR_ver}
  mkdir build
  cd build
  # build compilation option list

  COMPILATION_OPTIONS="-DCMAKE_INSTALL_PREFIX=${pkg_install_dir} -DUSE_OPENMP=ON -DBUILD_TESTING=NO -DWITH_EXAMPLES=NO"
  if [ "$MPI_MODE" == "no" ]; then
    COMPILATION_OPTIONS="${COMPILATION_OPTIONS} -DUSE_MPI=no"
  fi
  if [ "$ENABLE_CUDA" == "__TRUE__" ]; then
    COMPILATION_OPTIONS="${COMPILATION_OPTIONS} -DUSE_ACCEL=nvidia -DWITH_GPU=P100"
  fi
  if [ "$ENABLE_HIP" == "__TRUE__" ]; then
    COMPILATION_OPTIONS="${COMPILATION_OPTIONS} -DUSE_ACCEL=hip -DWITH_GPU=Mi250"
  fi
  cmake $COMPILATION_OPTIONS ..
  make -j $(get_nprocs) #> make.log 2>&1
  make install          #> install.log 2>&1
  cd ..
  write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
  DBCSR_CFLAGS="-I'${pkg_install_dir}/include'"
  DBCSR_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
fi
###        ;;
###    __SYSTEM__)
###        echo "==================== Finding DBCSR from system paths ===================="
###        check_lib -ldbcsr "dbcsr"
###        add_include_from_paths DBCSR_CFLAGS "dbcsr.h" $INCLUDE_PATHS
###        add_lib_from_paths DBCSR_LDFLAGS "dbcsr.*" $LIB_PATHS
###        ;;
###    __DONTUSE__)
###        report_error ${LINENO} "It is not possible to compile cp2k without dbcsr"
###        ;;
###    *)
echo "==================== Linking DBCSR to user paths ===================="
#####pkg_install_dir="$with_dbcsr"
check_dir "${pkg_install_dir}/lib"
check_dir "${pkg_install_dir}/include"
DBCSR_CFLAGS="-I'${pkg_install_dir}/include'"
DBCSR_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
##;;
##esac
##if [ "$with_dbcsr" != "__DONTUSE__" ]; then
DBCSR_LIBS="-ldbcsr"
##    if [ "$with_dbscr" != "__SYSTEM__" ]; then
cat << EOF > "${BUILDDIR}/setup_dbcsr"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path CMAKE_INSTALL_PREFIX "${pkg_install_dir}"
EOF
##        cat "${BUILDDIR}/setup_dbcsr" >> $SETUPFILE
##    fi
cat << EOF >> "${BUILDDIR}/setup_dbcsr"
export DBCSR_CFLAGS="${DBCSR_CFLAGS}"
export DBCSR_LDFLAGS="${DBCSR_LDFLAGS}"
export DBCSR_LIBS="${DBCSR_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_CUDA(-D__DBCSR_ACC -D__DBCSR|IF_HIP(-D__DBCSR_ACC -D__DBCSR|-D__DBCSR))"
export CP_CFLAGS="\${CP_CFLAGS} ${DBCSR_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${DBCSR_LDFLAGS}"
export CP_LIBS="${DBCSR_LIBS} \${CP_LIBS}"
export DBCSR_ROOT="$pkg_install_dir"
EOF

cat "${BUILDDIR}/setup_dbcsr" >> $SETUPFILE

#else
#    echo "DBCSR is a hard dependency for cp2k"
#    exit 1
#fi

load "${BUILDDIR}/setup_dbcsr"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "DBCSR"
