#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

json_fortran_ver=${json_fortran_ver:-7.0.0}
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_json_fortran=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_json_fortran" ] && rm -f "${BUILDDIR}/setup_json_fortran"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_json_fortran" in
    __INSTALL__)
        echo "==================== Installing json_fortran ===================="
        pkg_install_dir="${INSTALLDIR}/json_fortran-${json_fortran_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "json-fortran-${json_fortran_ver} is already installed, skipping it."
        else
            if [ -f json-fortran-${json_fortran_ver}.tar.gz ] ; then
                echo "json-fortran-${json_fortran_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://github.com/jacobwilliams/json-fortran/archive/${json_fortran_ver}.tar.gz \
                             -o json-fortran-${json_fortran_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d json-fortran-${json_fortran_ver} ] && rm -rf json-fortran-${json_fortran_ver}
            tar -xzf json-fortran-${json_fortran_ver}.tar.gz
            cd json-fortran-${json_fortran_ver}
            mkdir build
            cd build
            cmake -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" -DUSE_GNU_INSTALL_CONVENTION=on -DSKIP_DOC_GEN=true -DCMAKE_INSTALL_LIBDIR=lib .. > make.log 2>&1
            make -j $NPROCS >> make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ../..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
        fi

        JSON_CFLAGS="-I'${pkg_install_dir}/include'"
        JSON_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding json-fortran from system paths ===================="
        add_include_from_paths JSON_CFLAGS "json_file_module.mod" $INCLUDE_PATHS
        add_lib_from_paths JSON_LDFLAGS "libjsonfortran.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking json-fortran to user paths ===================="
        pkg_install_dir="$with_json"
        check_dir "$pkg_install_dir/lib"
        check_dir "$pkg_install_dir/lib64"
        check_dir "$pkg_install_dir/include"
        JSON_CFLAGS="-I'${pkg_install_dir}/include'"
        JSON_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_json_fortran" != "__DONTUSE__" ] ; then
    JSON_LIBS="-ljsonfortran"
    if [ "$with_json_fortran" != "__SYSTEM__" ] ; then
        cat << EOF > "${BUILDDIR}/setup_json_fortran"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
export JSON_CFLAGS="${JSON_CFLAGS}"
export JSON_LDFLAGS="${JSON_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} -D__JSON"
export CP_CFLAGS="\${CP_CFLAGS} ${JSON_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${JSON_LDFLAGS}"
export CP_LIBS="${JSON_LIBS} \${CP_LIBS}"
EOF
        cat "${BUILDDIR}/setup_json_fortran" >> $SETUPFILE
    fi
fi
cd "${ROOTDIR}"
