#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_scalapack=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_scalapack" ] && rm "${BUILDDIR}/setup_scalapack"

SCALAPACK_CFLAGS=''
SCALAPACK_LDFLAGS=''
SCALAPACK_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_scalapack" in
    __INSTALL__)
        echo "==================== Installing ScaLAPACK ===================="
        pkg_install_dir="${INSTALLDIR}/scalapack-${scalapack_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "scalapack-${scalapack_ver} is already installed, skipping it."
        else
            require_env MATH_LIBS
            if [ -f scalapack-${scalapack_ver}.tgz ] ; then
                echo "scalapack-${scalapack_ver}.tgz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/scalapack-${scalapack_ver}.tgz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d scalapack-${scalapack_ver} ] && rm -rf scalapack-${scalapack_ver}
            tar -xzf scalapack-${scalapack_ver}.tgz
            cd scalapack-${scalapack_ver}
            cat << EOF > SLmake.inc
CDEFS         = -DAdd_
FC            = ${MPIFC}
CC            = ${MPICC}
NOOPT         = ${FFLAGS} -O0 -fno-fast-math
FCFLAGS       = ${FFLAGS} ${MATH_CFLAGS}
CCFLAGS       = ${CFLAGS} ${MATH_CFLAGS}
FCLOADER      = \$(FC)
CCLOADER      = \$(CC)
FCLOADFLAGS   = \$(FCFLAGS) -Wl,--enable-new-dtags ${MATH_LDFLAGS}
CCLOADFLAGS   = \$(CCFLAGS) -Wl,--enable-new-dtags ${MATH_LDFLAGS}
ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = ranlib
SCALAPACKLIB  = libscalapack.a
BLASLIB       =
LAPACKLIB     = ${MATH_LIBS}
LIBS          = \$(LAPACKLIB) \$(BLASLIB)
EOF
            # scalapack build not parallel safe (update to the archive race)
            make -j 1 lib > make.log 2>&1
            # does not have make install, so install manually
            ! [ -d "${pkg_install_dir}/lib" ] && mkdir -p "${pkg_install_dir}/lib"
            cp libscalapack.a "${pkg_install_dir}/lib"
            cd ..
            touch "${install_lock_file}"
        fi
        SCALAPACK_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding ScaLAPACK from system paths ===================="
        check_lib -lscalapack "ScaLAPACK"
        add_lib_from_paths SCALAPACK_LDFLAGS "libscalapack.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking ScaLAPACK to user paths ===================="
        pkg_install_dir="$with_scalapack"
        check_dir "${pkg_install_dir}/lib"
        SCALAPACK_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_scalapack" != "__DONTUSE__" ] ; then
    SCALAPACK_LIBS="-lscalapack"
    if [ "$with_scalapack" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_scalapack"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
EOF
        cat "${BUILDDIR}/setup_scalapack" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_scalapack"
export SCALAPACK_LDFLAGS="${SCALAPACK_LDFLAGS}"
export SCALAPACK_LIBS="${SCALAPACK_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__SCALAPACK|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${SCALAPACK_LDFLAGS}|)"
export CP_LIBS="IF_MPI(-lscalapack|) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
