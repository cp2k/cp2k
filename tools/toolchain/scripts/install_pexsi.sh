#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

pexsi_ver="1.2.0"
pexsi_sha256="8bfad6ec6866c6a29e1cc87fb1c17a39809795e79ede98373c8ba9a3aaf820dd"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_pexsi" ] && rm "${BUILDDIR}/setup_pexsi"

PEXSI_CFLAGS=''
PEXSI_LDFLAGS=''
PEXSI_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_pexsi" in
    __INSTALL__)
        echo "==================== Installing PEXSI ===================="
        require_env SUPERLU_LDFLAGS
        require_env SUPERLU_LIBS
        require_env MATH_LIBS
        require_env MPI_LDFLAGS
        require_env MPI_LIBS
        pkg_install_dir="${INSTALLDIR}/pexsi-${pexsi_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "pexsi_dist-${pexsi_ver} is already installed, skipping it."
        else
            if [ -f pexsi_v${pexsi_ver}.tar.gz ] ; then
                echo "pexsi_v${pexsi_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} ${pexsi_sha256} \
                             https://www.cp2k.org/static/downloads/pexsi_v${pexsi_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d pexsi_v${pexsi_ver} ] && rm -rf pexsi_v${pexsi_ver}
            tar -xzf pexsi_v${pexsi_ver}.tar.gz
            cd pexsi_v${pexsi_ver}
            cat config/make.inc.linux.gnu | \
                sed -e "s|\(CC *=\).*|\1 ${MPICC}|g" \
                    -e "s|\(CXX *=\).*|\1 ${MPICXX}|g" \
                    -e "s|\(FC *=\).*|\1 ${MPIFC}|g" \
                    -e "s|\(LOADER *=\).*|\1 ${MPICXX}|g" \
                    -e "s|\(USE_SYMPACK *=\).*|\1 0|g" \
                    -e "s|\(PEXSI_DIR *=\).*|\1 ${PWD}|g" \
                    -e "s|\(PEXSI_BUILD_DIR *=\).*|\1 ${pkg_install_dir}|g" \
                    -e "s|\(CPP_LIB *=\).*|\1 -lstdc++ ${MPI_LDFLAGS} ${MPI_LIBS} |g" \
                    -e "s|\(LAPACK_LIB *=\).*|\1 ${MATH_LDFLAGS} $(resolve_string "${MATH_LIBS}")|g" \
                    -e "s|\(BLAS_LIB *=\).*|\1|g" \
                    -e "s|\(\bMETIS_LIB *=\).*|\1 ${METIS_LDFLAGS} ${METIS_LIBS}|g" \
                    -e "s|\(PARMETIS_LIB *=\).*|\1 |g" \
                    -e "s|\(DSUPERLU_LIB *=\).*|\1 ${SUPERLU_LDFLAGS} -lsuperlu_dist|g" \
                    -e "s|\(SCOTCH_LIB *=\).*|\1 ${SCOTCH_LDFLAGS} -lscotchmetis -lscotch -lscotcherr|g" \
                    -e "s|\(PTSCOTCH_LIB *=\).*|\1 ${SCOTCH_LDFLAGS} -lptscotchparmetis -lptscotch -lptscotcherr -lscotch|g" \
                    -e "s|\(DSUPERLU_INCLUDE *=\).*|\1 ${SUPERLU_CFLAGS}|g" \
                    -e "s|\(INCLUDES *=\).*|\1 ${MATH_CFLAGS} \${DSUPERLU_INCLUDE} \${PEXSI_INCLUDE}|g" \
                    -e "s|\(COMPILE_FLAG *=\).*|\1 ${CFLAGS}|g" \
                    -e "s|\(SUFFIX *=\).*|\1 ${OPENBLAS_ARCH}|g" \
                    -e "s|\(DSUPERLU_DIR *=\).*|\1|g" \
                    -e "s|\(METIS_DIR *=\).*|\1|g" \
                    -e "s|\(PARMETIS_DIR *=\).*|\1|g" \
                    -e "s|\(PTSCOTCH_DIR *=\).*|\1|g" \
                    -e "s|\(LAPACK_DIR *=\).*|\1|g" \
                    -e "s|\(BLAS_DIR *=\).*|\1|g" \
                    -e "s|-DCOREDUMPER||g" \
                    -e "s|\(COREDUMPER_LIB *=\).*||g" \
                    -e "s|\(GFORTRAN_LIB *=\).*|\1 -lgfortran|g" > make.inc

            # Seemingly this should've been moved together with the other fortran files.
            mv ./src/f_interface.f90 ./fortran/

            make finstall > make.log 2>&1 # still issues with parallel make (fortran_examples target)

            ln -sf "${pkg_install_dir}/lib/libpexsi_${OPENBLAS_ARCH}.a" \
                   "${pkg_install_dir}/lib/libpexsi.a"

            cp -r ./include/* ${pkg_install_dir}/include/  # bug: make install neglects most header files

            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi
        PEXSI_CFLAGS="-I'${pkg_install_dir}/include'"

        PEXSI_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding Pexsi_DIST from system paths ===================="
        check_lib -lpexsi "PEXSI"
        # add_include_from_paths PEXSI_CFLAGS "pexsi*" $INCLUDE_PATHS
        add_lib_from_paths PEXSI_LDFLAGS "libpexsi.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking Pexsi_Dist to user paths ===================="
        pkg_install_dir="$with_pexsi"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        ;;
esac
if [ "$with_pexsi" != "__DONTUSE__" ] ; then
    PEXSI_LIBS="-lpexsi"
    if [ "$with_pexsi" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_pexsi"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_pexsi" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_pexsi"
export PEXSI_CFLAGS="${PEXSI_CFLAGS}"
export PEXSI_LDFLAGS="${PEXSI_LDFLAGS}"
export PEXSI_LIBS="${PEXSI_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__LIBPEXSI|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(${PEXSI_CFLAGS}|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${PEXSI_LDFLAGS}|)"
export CP_LIBS="IF_MPI(${PEXSI_LIBS}|) \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_pexsi"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "pexsi"
