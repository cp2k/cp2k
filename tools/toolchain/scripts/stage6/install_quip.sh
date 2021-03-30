#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

quip_ver="b4336484fb65b0e73211a8f920ae4361c7c353fd"
quip_sha256="60fe54d60f5bcccd99abdccb6ca8d5d59c3c1c6997f95cee775318137743084e"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_quip" ] && rm "${BUILDDIR}/setup_quip"

if [ "${ENABLE_TSAN}" = "__TRUE__" ]; then
  report_warning "QUIP is not combatiable with thread sanitizer, not installing..."
  cat << EOF > setup_quip
with_quip=__DONTUSE__
EOF
  exit 0
fi

QUIP_CFLAGS=''
QUIP_LDFLAGS=''
QUIP_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_quip" in
  __INSTALL__)
    echo "==================== Installing QUIP ===================="
    require_env MATH_LIBS
    pkg_install_dir="${INSTALLDIR}/quip-${quip_ver}"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "quip_dist-${quip_ver} is already installed, skipping it."
    else
      if [ -f QUIP-${quip_ver}.tar.gz ]; then
        echo "QUIP-${quip_ver}.tar.gz is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} ${quip_sha256} \
          https://www.cp2k.org/static/downloads/QUIP-${quip_ver}.tar.gz
      fi
      [ -d QUIP-${quip_ver} ] && rm -rf QUIP-${quip_ver}
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf QUIP-${quip_ver}.tar.gz
      cd QUIP-${quip_ver}
      # translate OPENBLAS_ARCH
      case $OPENBLAS_ARCH in
        x86_64)
          quip_arch=x86_64
          ;;
        i386)
          quip_arch=x86_32
          ;;
        arm)
          quip_arch=x86_64
          ;;
        *)
          report_error $LINENO "arch $OPENBLAS_ARCH is currently unsupported"
          exit 1
          ;;
      esac
      # The ARCHER cd has a very annoying habbit of printing out
      # dir names to stdout for any target directories that are
      # more than one level deep, and one cannot seem to disable
      # it. This unfortunately messes up the installation script
      # for QUIP. So this hack will help to resolve the issue
      if [ "$ENABLE_CRAY" = "__TRUE__" ]; then
        sed -i \
          -e "s|\(cd build/.*\)|\1 >&- 2>&-|g" \
          bin/find_sizeof_fortran_t
      fi
      sed -i \
        -e "s|\(F77 *=\).*|\1 ${FC}|g" \
        -e "s|\(F90 *=\).*|\1 ${FC}|g" \
        -e "s|\(F95 *=\).*|\1 ${FC}|g" \
        -e "s|\(CC *=\).*|\1 ${CC}|g" \
        -e "s|\(CPLUSPLUS *=\).*|\1 ${CXX}|g" \
        -e "s|\(LINKER *=\).*|\1 ${FC}|g" \
        -e "s|\(FPP *=\).*|\1 ${FC} -E -x f95-cpp-input|g" \
        -e "s|\(QUIPPY_FCOMPILER *=\).*|\1 ${FC}|g" \
        -e "s|\(QUIPPY_CPP *=\).*|\1 ${FC} -E -x f95-cpp-input|g" \
        arch/Makefile.linux_${quip_arch}_gfortran

      # workaround for compilation with GCC-10, until properly fixed:
      #   https://github.com/libAtoms/QUIP/issues/209
      ("${FC}" --version | grep -Eq 'GNU.+\s10\.') && compat_flag="-fallow-argument-mismatch" || compat_flag=""

      # enable debug symbols
      echo "F95FLAGS       += -g ${compat_flag}" >> arch/Makefile.linux_${quip_arch}_gfortran
      echo "F77FLAGS       += -g ${compat_flag}" >> arch/Makefile.linux_${quip_arch}_gfortran
      echo "CFLAGS         += -g" >> arch/Makefile.linux_${quip_arch}_gfortran
      echo "CPLUSPLUSFLAGS += -g" >> arch/Makefile.linux_${quip_arch}_gfortran
      export QUIP_ARCH=linux_${quip_arch}_gfortran
      # hit enter a few times to accept defaults
      echo -e "${MATH_LDFLAGS} $(resolve_string "${MATH_LIBS}") \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" | make config > configure.log
      # make -j does not work :-(
      make > make.log 2>&1
      ! [ -d "${pkg_install_dir}/include" ] && mkdir -p "${pkg_install_dir}/include"
      ! [ -d "${pkg_install_dir}/lib" ] && mkdir -p "${pkg_install_dir}/lib"
      cp build/linux_x86_64_gfortran/quip_unified_wrapper_module.mod \
        "${pkg_install_dir}/include/"
      cp build/linux_x86_64_gfortran/*.a \
        "${pkg_install_dir}/lib/"
      cp src/fox/objs.linux_${quip_arch}_gfortran/lib/*.a \
        "${pkg_install_dir}/lib/"
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    QUIP_CFLAGS="-I'${pkg_install_dir}/include'"
    QUIP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding Quip_DIST from system paths ===================="
    check_lib -lquip_core "QUIP"
    check_lib -latoms "QUIP"
    check_lib -lFoX_sax "QUIP"
    check_lib -lFoX_common "QUIP"
    check_lib -lFoX_utils "QUIP"
    check_lib -lFoX_fsys "QUIP"
    add_include_from_paths QUIP_CFLAGS "quip*" $INCLUDE_PATHS
    add_lib_from_paths QUIP_LDFLAGS "libquip_core*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Quip_Dist to user paths ===================="
    pkg_install_dir="$with_quip"
    check_dir "${pkg_install_dir}/lib"
    check_dir "${pkg_install_dir}/include"
    QUIP_CFLAGS="-I'${pkg_install_dir}/include'"
    QUIP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_quip" != "__DONTUSE__" ]; then
  QUIP_LIBS="-lquip_core -latoms -lFoX_sax -lFoX_common -lFoX_utils -lFoX_fsys"
  if [ "$with_quip" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_quip"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
    cat "${BUILDDIR}/setup_quip" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_quip"
export QUIP_CFLAGS="${QUIP_CFLAGS}"
export QUIP_LDFLAGS="${QUIP_LDFLAGS}"
export QUIP_LIBS="${QUIP_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__QUIP"
export CP_CFLAGS="\${CP_CFLAGS} ${QUIP_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${QUIP_LDFLAGS}"
export CP_LIBS="${QUIP_LIBS} \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_quip"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "quip"
