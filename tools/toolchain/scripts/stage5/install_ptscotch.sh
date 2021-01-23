#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

scotch_ver="6.0.0"
scotch_sha256="e57e16c965bab68c1b03389005ecd8a03745ba20fd9c23081c0bb2336972d879"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_ptscotch" ] && rm "${BUILDDIR}/setup_ptscotch"

SCOTCH_CFLAGS=''
SCOTCH_LDFLAGS=''
SCOTCH_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_ptscotch" in
    __INSTALL__)
        echo "==================== Installing PT-Scotch ===================="
        pkg_install_dir="${INSTALLDIR}/scotch-${scotch_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "scotch-${scotch_ver} is already installed, skipping it."
        else
            if [ -f scotch_${scotch_ver}.tar.gz ] ; then
                echo "scotch_${scotch_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${scotch_sha256} \
                             https://www.cp2k.org/static/downloads/scotch_${scotch_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d scotch_${scotch_ver} ] && rm -rf scotch_${scotch_ver}
            tar -xzf scotch_${scotch_ver}.tar.gz
            cd scotch_${scotch_ver}/src
            cat Make.inc/Makefile.inc.x86-64_pc_linux2 | \
                sed -e "s|\(^CCS\).*|\1 = $MPICC|g" \
                    -e "s|\(^CCP\).*|\1 = $MPICC|g" \
                    -e "s|\(^CCD\).*|\1 = $MPICC|g" \
                    -e "s|\(^CFLAGS\).*|\1 = $CFLAGS -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict -DIDXSIZE64 ${MPI_CFLAGS}|g" \
                    > Makefile.inc
            make scotch -j $NPROCS > make.log 2>&1
            make ptscotch -j $NROCS > make.log 2>&1
            # PT-scotch make install is buggy in that it cannot create
            # intermediate directories
            ! [ -d "${pkg_install_dir}" ] && mkdir -p "${pkg_install_dir}"
            make install prefix=${pkg_install_dir} > install.log 2>&1
            cd ..

            # PEXSI also needs parmetis.h
            cp ./include/parmetis.h "${pkg_install_dir}/include/"
            sed -i "s|SCOTCH_Num|int|g" "${pkg_install_dir}/include/parmetis.h"

            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage5/$(basename ${SCRIPT_NAME})"
        fi
        SCOTCH_CFLAGS="-I'${pkg_install_dir}/include'"
        SCOTCH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding PT-Scotch from system paths ===================="
        check_lib -lptscotchparmetis "PT-Scotch"
        check_lib -lptscotch "PT-Scotch"
        check_lib -lptscotcherr "PT-Scotch"
        check_lib -lscotchmetis "PT-Scotch"
        check_lib -lscotch "PT-Scotch"
        check_lib -lscotcherr "PT-Scotch"
        check_lib -lptscotchparmetis "PT-Scotch"
        add_include_from_paths SCOTCH_CFLAGS "ptscotch.h" $INCLUDE_PATHS
        add_lib_from_paths SCOTCH_LDFLAGS "libptscotch.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking PT-Scotch to user paths ===================="
        pkg_install_dir="$with_ptscotch"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        SCOTCH_CFLAGS="-I'${pkg_install_dir}/include'"
        SCOTCH_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_ptscotch" != "__DONTUSE__" ] ; then
    SCOTCH_LIBS="-lptscotchparmetis -lptscotch -lptscotcherr -lscotchmetis -lscotch -lscotcherr"
    if [ "$with_ptscotch" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_ptscotch"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_ptscotch" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_ptscotch"
export SCOTCH_CFLAGS="${SCOTCH_CFLAGS}"
export SCOTCH_LDFLAGS="${SCOTCH_LDFLAGS}"
export SCOTCH_LIBS="${SCOTCH_LIBS}"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${SCOTCH_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${SCOTCH_LDFLAGS}|)"
export CP_LIBS="IF_MPI(${SCOTCH_LIBS}|) \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_ptscotch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "ptscotch"
