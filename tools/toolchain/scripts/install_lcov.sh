#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_lcov=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_lcov" ] && rm "${BUILDDIR}/setup_lcov"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_lcov" in
    __INSTALL__)
        echo "==================== Installing Lcov ===================="
        pkg_install_dir="${INSTALLDIR}/lcov-${lcov_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "lcov-${lcov_ver} is already installed, skipping it."
        else
            if [ -f lcov-${lcov_ver}.tar.gz ] ; then
                echo "lcov-${lcov_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/lcov-${lcov_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d lcov-${lcov_ver} ] && rm -rf lcov-${lcov_ver}
            tar -xzf lcov-${lcov_ver}.tar.gz
            cd lcov-${lcov_ver}
            # note.... this installs in ${INSTALLDIR}/usr/bin
            make PREFIX="${pkg_install_dir}" install > make.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        ;;
    __SYSTEM__)
        echo "==================== Finding Lcov from system paths ===================="
        check_command lcov "lcov"
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking Lcov to user paths ===================="
        pkg_install_dir="$with_lcov"
        check_dir "${with_lcov}/usr/bin"
        ;;
esac
if [ "$with_lcov" != "__DONTUSE__" ] ; then
    if [ "$with_lcov" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_lcov"
prepend_path PATH "$pkg_install_dir/usr/bin"
EOF
        cat "${BUILDDIR}/setup_lcov" >> ${SETUPFILE}
    fi
fi
cd "${ROOTDIR}"
