#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh

with_valgrind=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_valgrind" ] && rm "${BUILDDIR}/setup_valgrind"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_valgrind" in
    __INSTALL__)
        echo "==================== Installing Valgrind ===================="
        pkg_install_dir="${INSTALLDIR}/valgrind-${valgrind_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "valgrind-${valgrind_ver} is already installed, skipping it."
        else
            if [ -f valgrind-${valgrind_ver}.tar.bz2 ] ; then
                echo "valgrind-${valgrind_ver}.tar.bz2 is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/valgrind-${valgrind_ver}.tar.bz2
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xjf valgrind-${valgrind_ver}.tar.bz2
            cd valgrind-${valgrind_ver}
            ./configure --prefix="${pkg_install_dir}" > config.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            touch "${install_lock_file}"
        fi
        ;;
    __SYSTEM__)
        echo "==================== Finding Valgrind from system paths ===================="
        check_command valgrind "valgrind"
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking Valgrind to user paths ===================="
        pkg_install_dir="$with_valgrind"
        check_dir "${with_valgrind}/bin"
        check_dir "${with_valgrind}/lib"
        check_dir "${with_valgrind}/include"
        ;;
esac
if [ "$with_valgrind" != "__DONTUSE__" ] ; then
    if [ "$with_valgrind" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_valgrind"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path PATH "$pkg_install_dir/lib"
prepend_path PATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_valgrind" >> ${SETUPFILE}
    fi
fi
cd "${ROOTDIR}"

# ----------------------------------------------------------------------
# Suppress reporting of known leaks
# ----------------------------------------------------------------------
cat <<EOF > ${INSTALLDIR}/valgrind.supp
{
   BuggySUPERLU
   Memcheck:Cond
   ...
   fun:SymbolicFactorize
}
{
   BuggyMPICH32
   Memcheck:Cond
   ...
   fun:MPIR_Process_status
}
{
   BuggyLD
   Memcheck:Cond
   ...
   fun:expand_dynamic_string_token
}
EOF
# also need to give links to the .supp file in setup file
cat <<EOF >> ${SETUPFILE}
export VALGRIND_OPTIONS="--suppressions=${INSTALLDIR}/valgrind.supp --max-stackframe=2168152 --error-exitcode=42"
EOF
