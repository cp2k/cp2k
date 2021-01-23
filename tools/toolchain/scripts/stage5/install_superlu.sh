#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

superlu_ver="6.1.0"  # Newer versions don't work with PEXSI 1.2.0.
superlu_sha256="92c6d1424dd830ee2d1e7396a418a5f6645160aea8472e558c4e4bfe006593c4"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_superlu" ] && rm "${BUILDDIR}/setup_superlu"

SUPERLU_CFLAGS=''
SUPERLU_LDFLAGS=''
SUPERLU_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_superlu" in
    __INSTALL__)
        echo "==================== Installing SuperLU_DIST ===================="
        require_env MATH_LIBS
        pkg_install_dir="${INSTALLDIR}/superlu_dist-${superlu_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "superlu_dist-${superlu_ver} is already installed, skipping it."
        else
            if [ -f superlu_dist_${superlu_ver}.tar.gz ] ; then
                echo "superlu_dist_${superlu_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${superlu_sha256} \
                             https://www.cp2k.org/static/downloads/superlu_dist_${superlu_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d SuperLU_DIST_${superlu_ver} ] && rm -rf SuperLU_DIST_${superlu_ver}
            tar -xzf superlu_dist_${superlu_ver}.tar.gz
            cd superlu_dist-${superlu_ver}
            cd build
            # Explicitly set LIBDIR to "lib", otherwise it sometimes defaults to "lib64".
            cmake -DTPL_ENABLE_PARMETISLIB=FALSE \
                  -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} \
                  -DCMAKE_INSTALL_LIBDIR=${pkg_install_dir}/lib \
                  .. > cmake.log 2>&1
            make -j $NPROCS > make.log 2>&1
            make install > install.log 2>&1
            cd ..

            # PEXSI needs some more headers.
            cp SRC/*.h "${pkg_install_dir}/include"

            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage5/$(basename ${SCRIPT_NAME})"
        fi
        SUPERLU_CFLAGS="-I'${pkg_install_dir}/include'"
        SUPERLU_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding SuperLU_DIST from system paths ===================="
        check_lib -lsuperlu_dist "SuperLU_DIST"
        add_include_from_paths SUPERLU_CFLAGS "superlu*" $INCLUDE_PATHS
        add_lib_from_paths SUPERLU_LDFLAGS "libsuperlu*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking Superlu_Dist to user paths ===================="
        pkg_install_dir="$with_superlu"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        SUPERLU_CFLAGS="-I'${pkg_install_dir}/include'"
        SUPERLU_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_superlu" != "__DONTUSE__" ] ; then
    SUPERLU_LIBS="-lsuperlu_dist"
    if [ "$with_superlu" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_superlu"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_superlu" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_superlu"
export SUPERLU_CFLAGS="${SUPERLU_CFLAGS}"
export SUPERLU_LDFLAGS="${SUPERLU_LDFLAGS}"
export SUPERLU_LIBS="${SUPERLU_LIBS}"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${SUPERLU_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${SUPERLU_LDFLAGS}|)"
export CP_LIBS="IF_MPI(${SUPERLU_LIBS}|) \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_superlu"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "superlu"
