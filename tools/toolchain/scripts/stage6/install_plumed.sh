#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

plumed_ver="2.6.2"
plumed_pkg="plumed-${plumed_ver}.tgz"
plumed_sha256="1ab3153db2010406852b30201ed94112e25eca4c4c8c4b41a29c22a7a3303f96"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_plumed" ] && rm "${BUILDDIR}/setup_plumed"

PLUMED_LDFLAGS=''
PLUMED_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_plumed" in
  __INSTALL__)
    echo "==================== Installing PLUMED ===================="
    pkg_install_dir="${INSTALLDIR}/plumed-${plumed_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "plumed-${plumed_ver} is already installed, skipping it."
    else
      if [ -f ${plumed_pkg} ]; then
        echo "${plumed_pkg} is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${plumed_sha256} \
          "https://www.cp2k.org/static/downloads/${plumed_pkg}"
      fi

      [ -d plumed-${plumed_ver} ] && rm -rf plumed-${plumed_ver}
      tar -xzf ${plumed_pkg}

      echo "Installing from scratch into ${pkg_install_dir}"
      cd plumed-${plumed_ver}
      # disable generating debugging infos for now to work around an issue in gcc-10.2:
      # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96354
      # note: some MPI wrappers carry a -g forward, thus stripping is not enough

      libs=""
      [ -n "${MKL_LIBS}" ] && libs+="${MKL_LIBS}"

      ./configure \
        CXX="${MPICXX}" \
        CXXFLAGS="${CXXFLAGS//-g/-g0} ${GSL_CFLAGS}" \
        LDFLAGS="${LDFLAGS} ${GSL_LDFLAGS}" \
        LIBS="${libs}" \
        --prefix=${pkg_install_dir} \
        --libdir="${pkg_install_dir}/lib" > configure.log 2>&1
      make VERBOSE=1 -j $NPROCS > make.log 2>&1
      make install > install.log 2>&1
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    PLUMED_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding PLUMED from system paths ===================="
    check_lib -lplumed "PLUMED"
    add_lib_from_paths PLUMED_LDFLAGS "libplumed*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking PLUMED to user paths ===================="
    pkg_install_dir="$with_plumed"
    check_dir "${pkg_install_dir}/lib"
    PLUMED_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac

if [ "$with_plumed" != "__DONTUSE__" ]; then
  PLUMED_LIBS='-lplumed -ldl -lstdc++ -lz -ldl'
  if [ "$with_plumed" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_plumed"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
EOF
    cat "${BUILDDIR}/setup_plumed" >> $SETUPFILE
  fi

  cat << EOF >> "${BUILDDIR}/setup_plumed"
export PLUMED_LDFLAGS="${PLUMED_LDFLAGS}"
export PLUMED_LIBS="${PLUMED_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__PLUMED2"
export CP_LDFLAGS="\${CP_LDFLAGS} ${PLUMED_LDFLAGS}"
export CP_LIBS="${PLUMED_LIBS} \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_plumed"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "plumed"
