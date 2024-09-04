#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# Installing QUIP without GAP because its ASL licence is not GPL-compatible.
# See also https://github.com/libAtoms/QUIP/issues/481

quip_ver="0.9.10"
quip_sha256="c03505779634459ea0ba3f7ddc120ac17f0546d44dc9b5096f008f1c3c6620ef"

fox_ver="b5b69ef9a46837bd944ba5c9bc1cf9d00a6198a7"
fox_sha256="a87dd7faf80612a0df94dc272474f37689c6213c6ac4705fb637644409c5cd4e"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_quip" ] && rm "${BUILDDIR}/setup_quip"

if [ "${ENABLE_TSAN}" = "__TRUE__" ]; then
  report_warning "QUIP is not compatible with the thread sanitizer. The QUIP package will not be installed."
  cat << EOF > ${BUILDDIR}/setup_quip
with_quip="__DONTUSE__"
EOF
  exit 0
fi

if [ "${with_quip}" != "__DONTUSE__" ] && [ "${with_intel}" != "__DONTUSE__" ]; then
  report_warning "A QUIP installation using the Intel compiler is currently not supported. The QUIP package will not be installed."
  exit 0
fi

QUIP_CFLAGS=""
QUIP_LDFLAGS=""
QUIP_LIBS=""
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_quip}" in
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
        download_pkg_from_cp2k_org "${quip_sha256}" "QUIP-${quip_ver}.tar.gz"
      fi
      if [ -f fox-${fox_ver}.tar.gz ]; then
        echo "fox-${fox_ver}.tar.gz is found"
      else
        download_pkg_from_cp2k_org "${fox_sha256}" "fox-${fox_ver}.tar.gz"
      fi
      [ -d QUIP-${quip_ver} ] && rm -rf QUIP-${quip_ver}
      [ -d fox-${fox_ver} ] && rm -rf fox-${fox_ver}
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf QUIP-${quip_ver}.tar.gz
      tar -xzf fox-${fox_ver}.tar.gz
      cd QUIP-${quip_ver}
      rmdir ./src/fox
      ln -s ../../fox-${fox_ver} ./src/fox
      # translate OPENBLAS_ARCH
      case $OPENBLAS_ARCH in
        x86_64)
          quip_arch="x86_64"
          ;;
        i386)
          quip_arch="x86_32"
          ;;
        arm*)
          quip_arch="x86_64"
          ;;
        *)
          report_error ${LINENO} "arch $OPENBLAS_ARCH is currently not supported."
          exit 1
          ;;
      esac
      # The ARCHER cd has a very annoying habit of printing out
      # dir names to stdout for any target directories that are
      # more than one level deep, and one cannot seem to disable
      # it. This unfortunately messes up the installation script
      # for QUIP. So this hack will help to resolve the issue
      if [ "${ENABLE_CRAY}" = "__TRUE__" ]; then
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
      if ("${FC}" --version | grep -q 'GNU'); then
        compat_flag=$(allowed_gfortran_flags "-fallow-argument-mismatch")
      fi

      # enable debug symbols
      echo "F95FLAGS       += -g ${compat_flag}" >> arch/Makefile.linux_${quip_arch}_gfortran
      echo "F77FLAGS       += -g ${compat_flag}" >> arch/Makefile.linux_${quip_arch}_gfortran
      echo "CFLAGS         += -g -fpermissive" >> arch/Makefile.linux_${quip_arch}_gfortran
      echo "CPLUSPLUSFLAGS += -g -fpermissive" >> arch/Makefile.linux_${quip_arch}_gfortran
      # Makefile.linux_${quip_arch}_gfortran_openmp includes Makefile.linux_${quip_arch}_gfortran
      export QUIP_ARCH=linux_${quip_arch}_gfortran_openmp
      # hit enter a few times to accept defaults
      echo -e "${MATH_LDFLAGS} $(resolve_string "${MATH_LIBS}") \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" | make config > configure.log
      # make -j does not work :-(
      make > make.log 2>&1 || tail -n ${LOG_LINES} make.log
      ! [ -d "${pkg_install_dir}/include" ] && mkdir -p "${pkg_install_dir}/include"
      ! [ -d "${pkg_install_dir}/lib" ] && mkdir -p "${pkg_install_dir}/lib"
      cp build/${QUIP_ARCH}/quip_unified_wrapper_module.mod \
        "${pkg_install_dir}/include/"
      cp build/${QUIP_ARCH}/*.a \
        "${pkg_install_dir}/lib/"
      cp src/fox/objs.${QUIP_ARCH}/lib/*.a \
        "${pkg_install_dir}/lib/"
      cd ..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    QUIP_CFLAGS="-I'${pkg_install_dir}/include'"
    QUIP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
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
    QUIP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath,'${pkg_install_dir}/lib'"
    ;;
esac
if [ "${with_quip}" != "__DONTUSE__" ]; then
  QUIP_LIBS="-lquip_core -latoms -lFoX_sax -lFoX_common -lFoX_utils -lFoX_fsys"
  if [ "${with_quip}" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_quip"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
prepend_path PKG_CONFIG_PATH "$pkg_install_dir/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "$pkg_install_dir"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_quip"
export QUIP_CFLAGS="${QUIP_CFLAGS}"
export QUIP_LDFLAGS="${QUIP_LDFLAGS}"
export QUIP_LIBS="${QUIP_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__QUIP"
export CP_CFLAGS="\${CP_CFLAGS} ${QUIP_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${QUIP_LDFLAGS}"
export CP_LIBS="${QUIP_LIBS} \${CP_LIBS}"
export QUIP_ROOT="${pkg_install_dir}"
EOF
  cat "${BUILDDIR}/setup_quip" >> ${SETUPFILE}
fi

load "${BUILDDIR}/setup_quip"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "quip"
