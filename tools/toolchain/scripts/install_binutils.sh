#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_binutils=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_binutils" ] && rm "${BUILDDIR}/setup_binutils"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_binutils" in
    __INSTALL__)
        echo "==================== Installing binutils ===================="
        pkg_install_dir="${INSTALLDIR}/binutils-${binutils_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "binutils-${binutils_ver} is already installed, skipping it."
        else
            if [ -f binutils-${binutils_ver}.tar.gz ] ; then
                echo "binutils-${binutils_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://ftp.gnu.org/gnu/binutils/binutils-${binutils_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xzf binutils-${binutils_ver}.tar.gz
            cd binutils-${binutils_ver}
            ./configure --prefix="${pkg_install_dir}" \
                        --enable-gold --enable-plugins > configure.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        ;;
    __SYSTEM__)
        echo "==================== Finding binutils from system paths ===================="
        check_command ar "gnu binutils"
        check_command ld "gnu binutils"
        check_command ranlib "gnu binutils"
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking binutils to user paths ===================="
        pkg_install_dir="$with_binutils"
        check_dir "$pkg_install_dir/bin"
        check_dir "$pkg_install_dir/lib"
        check_dir "$pkg_install_dir/lib64"
        check_dir "$pkg_install_dir/include"
        ;;
esac
if [ "$with_binutils" != "__DONTUSE__" ] ; then
    if [ "$with_binutils" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_binutils"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib64"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib64"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib64"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_binutils" >> $SETUPFILE
    fi
fi
cd "${ROOTDIR}"
