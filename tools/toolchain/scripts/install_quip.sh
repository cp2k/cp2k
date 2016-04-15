#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_quip=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_quip" ] && rm "${BUILDDIR}/setup_quip"

if [ "${ENABLE_TSAN}" = "__TRUE__" ] ; then
    report_warning "QUIP is not combatiable with thread sanitizer, not installing..."
cat <<EOF > setup_quip
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
        if [ -f "${install_lock_file}" ] ; then
            echo "quip_dist-${quip_ver} is already installed, skipping it."
        else
            if [ -f QUIP-${quip_ver}.zip ] ; then
                echo "QUIP-${quip_ver}.zip is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/QUIP-${quip_ver}.zip
            fi
            [ -d QUIP-${quip_ver} ] && rm -rf QUIP-${quip_ver}
            echo "Installing from scratch into ${pkg_install_dir}"
            unzip -q -o QUIP-${quip_ver}.zip
            cd QUIP-${quip_ver}
            # translate OPENBLAS_ARCH
            case $OPENBLAS_ARCH in
                x86_64)
                    quip_arch=x86_64
                    ;;
                i386)
                    quip_arch=x86_32
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
            if [ "$ENABLE_CRAY" = "__TRUE__" ] ; then
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
            # enable debug symbols
            echo "F95FLAGS       += -g" >> arch/Makefile.linux_${quip_arch}_gfortran
            echo "F77FLAGS       += -g" >> arch/Makefile.linux_${quip_arch}_gfortran
            echo "CFLAGS         += -g" >> arch/Makefile.linux_${quip_arch}_gfortran
            echo "CPLUSPLUSFLAGS += -g" >> arch/Makefile.linux_${quip_arch}_gfortran
            export QUIP_ARCH=linux_${quip_arch}_gfortran
            # hit enter a few times to accept defaults
            echo -e "${MATH_LDFLAGS} ${MATH_LIBS} \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" | make config > configure.log
            # make -j does not work :-(
            make > make.log 2>&1
            ! [ -d "${pkg_install_dir}/include" ] && mkdir -p "${pkg_install_dir}/include"
            ! [ -d "${pkg_install_dir}/lib" ] && mkdir -p "${pkg_install_dir}/lib"
            cp build/linux_x86_64_gfortran/quip_unified_wrapper_module.mod \
               "${pkg_install_dir}/include/"
            cp build/linux_x86_64_gfortran/*.a \
               "${pkg_install_dir}/lib/"
            cp src/FoX-4.0.3/objs.linux_${quip_arch}_gfortran/lib/*.a \
               "${pkg_install_dir}/lib/"
            cd ..
            touch "${install_lock_file}"
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
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking Quip_Dist to user paths ===================="
        pkg_install_dir="$with_quip"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        QUIP_CFLAGS="-I'${pkg_install_dir}/include'"
        QUIP_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_quip" != "__DONTUSE__" ] ; then
    QUIP_LIBS="-lquip_core -latoms -lFoX_sax -lFoX_common -lFoX_utils -lFoX_fsys"
    if [ "$with_quip" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_quip"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_quip" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_quip"
export QUIP_CFLAGS="${QUIP_CFLAGS}"
export QUIP_LDFLAGS="${QUIP_LDFLAGS}"
export QUIP_LIBS="${QUIP_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__QUIP"
export CP_CFLAGS="\${CP_CFLAGS} ${QUIP_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${QUIP_LDFLAGS}"
export CP_LIBS="${QUIP_LIBS} \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
