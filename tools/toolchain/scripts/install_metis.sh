#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_metis=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_metis" ] && rm "${BUILDDIR}/setup_metis"

METIS_CFLAGS=''
METIS_LDFLAGS=''
METIS_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_metis" in
    __INSTALL__)
        echo "==================== Installing METIS ===================="
        require_env PARMETIS_INSTALL_MODE
        echo "METIS is installed together with ParMETIS"
        if [ "$PARMETIS_INSTALL_MODE" = __INSTALL__ ] ; then
            METIS_CFLAGS="$PARMETIS_CFLAGS"
            METIS_LDFLAGS="$PARMETIS_LDFLAGS"
        else
            report error $LINENO "Use option --with-parmetis=install to install METIS"
            exit 1
        fi
        ;;
    __SYSTEM__)
        echo "==================== Finding METIS from system paths ===================="
        check_lib -lmetis "METIS"
        add_include_from_paths METIS_CFLAGS "metis.h" $INCLUDE_PATHS
        add_lib_from_paths METIS_LDFLAGS "libmetis.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking METIS to user paths ===================="
        pkg_install_dir="$with_metis"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        METIS_CFLAGS="-I'${pkg_install_dir}/include'"
        METIS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_metis" != "__DONTUSE__" ] ; then
    METIS_LIBS="-lmetis"
    # we assume if in __INSTALL__ mode, then METIS should have been
    # installed from ParMETIS, and this means we don't need to setup
    # the system paths for METIS
    if [ "$with_metis" != "__SYSTEM__" ] && \
       [ "$with_metis" != "__INSTALL__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_metis"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_metis" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_metis"
export METIS_CFLAGS="${METIS_CFLAGS}"
export METIS_LDFLAGS="${METIS_LDFLAGS}"
export METIS_LIBS="${METIS_LIBS}"
EOF
fi
cd "${ROOTDIR}"
