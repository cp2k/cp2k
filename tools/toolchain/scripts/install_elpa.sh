#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/package_versions.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_elpa=${1:-__INSTALL__}

[ -f "${BUILDDIR}/setup_elpa" ] && rm "${BUILDDIR}/setup_elpa"

ELPA_CFLAGS=''
ELPA_LDFLAGS=''
ELPA_LIBS=''
ELPA_CFLAGS_OMP=''
ELPA_LIBS_OMP=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

# elpa only works with MPI switched on
if [ $MPI_MODE = no ] ; then
    report_warning $LINENO "MPI is disabled, skipping elpa installation"
cat <<EOF > "${BUILDDIR}/setup_elpa"
with_elpa="__FALSE__"
EOF
    exit 0
fi

case "$with_elpa" in
    __INSTALL__)
        echo "==================== Installing ELPA ===================="
        pkg_install_dir="${INSTALLDIR}/elpa-${elpa_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "elpa-${elpa_ver} is already installed, skipping it."
        else
            require_env MATH_LIBS
            if [ -f elpa-${elpa_ver}.tar.gz ] ; then
                echo "elpa-${elpa_ver}.tar.gz is found"
            else
                download_pkg ${DOWNLOADER_FLAGS} \
                             https://www.cp2k.org/static/downloads/elpa-${elpa_ver}.tar.gz
            fi
            [ -d elpa-${elpa_ver} ] && rm -rf elpa-${elpa_ver}
            echo "Installing from scratch into ${pkg_install_dir}"
            tar -xzf elpa-${elpa_ver}.tar.gz

            # fix wrong dependency order (at least in elpa 2016.05.003)
            sed -i "s/build_lib = libelpa@SUFFIX@.la libelpatest@SUFFIX@.la/build_lib = libelpatest@SUFFIX@.la libelpa@SUFFIX@.la/g" elpa-${elpa_ver}/Makefile.in
            sed -i "s/build_lib = libelpa@SUFFIX@.la libelpatest@SUFFIX@.la/build_lib = libelpatest@SUFFIX@.la libelpa@SUFFIX@.la/g" elpa-${elpa_ver}/Makefile.am

            # need both flavors ?

            # elpa expect FC to be an mpi fortran compiler that is happy
            # with long lines, and that a bunch of libs can be found
            cd elpa-${elpa_ver}
            # shared libraries cannot be built if linked with
            # reflapack (the case when valgrind is enabled)
            if [ $ENABLE_VALGRIND = "__TRUE__" ] ; then
                shared_flag=no
            else
                shared_flag=yes
            fi
            # specific settings needed on CRAY Linux Environment
            if [ $ENABLE_CRAY = "__TRUE__" ] ; then
                # extra LDFLAGS needed
                cray_ldflags="-dynamic"
            fi
            # non-threaded version
            ./configure  --prefix=${pkg_install_dir} \
                         --libdir="${pkg_install_dir}/lib" \
                         --enable-openmp=no \
                         --enable-shared=$shared_flag \
                         --enable-static=yes \
                         FC=${MPIFC} \
                         CC=${MPICC} \
                         CXX=${MPICXX} \
                         FCFLAGS="${FCFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} -ffree-line-length-none" \
                         CFLAGS="${CFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS}" \
                         CXXFLAGS="${CXXFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS}" \
                         LDFLAGS="-Wl,--enable-new-dtags ${MATH_LDFLAGS} ${SCALAPACK_LDFLAGS} ${cray_ldflags}" \
                         LIBS="${SCALAPACK_LIBS} ${MATH_LIBS}" \
                         > configure.log 2>&1
            make -j $NPROCS >  make.log 2>&1
            make install > install.log 2>&1
            # threaded version
            if [ "$ENABLE_OMP" = "__TRUE__" ] ; then
                make -j $NPROCS clean > clean.log 2>&1
                ./configure  --prefix=${pkg_install_dir} \
                             --enable-openmp=yes \
                             --enable-shared=$shared_flag \
                             --enable-static=yes \
                             FC=${MPIFC} \
                             CC=${MPICC} \
                             CXX=${MPICXX} \
                             FCFLAGS="${FCFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS} -ffree-line-length-none" \
                             CFLAGS="${CFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS}" \
                             CXXFLAGS="${CXXFLAGS} ${MATH_CFLAGS} ${SCALAPACK_CFLAGS}" \
                             LDFLAGS="-Wl,--enable-new-dtags ${MATH_LDFLAGS} ${SCALAPACK_LDFLAGS} ${cray_ldflags}" \
                             LIBS="${SCALAPACK_LIBS} ${MATH_LIBS}" \
                             > configure.log 2>&1
                make -j $NPROCS >  make.log 2>&1
                make install > install.log 2>&1
            fi
            cd ..
            touch "${install_lock_file}"
        fi
        ELPA_CFLAGS="-I'${pkg_install_dir}/include/elpa-${elpa_ver}/modules'"
        ELPA_CFLAGS_OMP="-I'${pkg_install_dir}/include/elpa_openmp-${elpa_ver}/modules'"
        ELPA_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
    __SYSTEM__)
        echo "==================== Finding ELPA from system paths ===================="
        check_lib -lelpa "ELPA"
        [ "$ENABLE_OMP" = "__TRUE__" ] && check_lib -lelpa_openmp "ELPA threaded version"
        # get the include paths
        elpa_include="$(find_in_paths "elpa-*" $INCLUDE_PATHS)"
        if [ "$elpa_include" != "__FALSE__" ] ; then
            echo "ELPA include directory is found to be $elpa_include/modules"
            elpa_include="$elpa_include/modules"
            ELPA_CFLAGS="-I'$elpa_include'"
        else
            echo "Cannot find elpa-* from paths $INCLUDE_PATHS"
            exit 1
        fi
        if [ "$ENABLE_OMP" = "__TRUE__" ] ; then
            elpa_include_omp="$(find_in_paths "elpa_openmp-*" $INCLUDE_PATHS)"
            if [ "$elpa_include_omp" != "__FALSE__" ] ; then
                echo "ELPA include directory threaded version is found to be $elpa_include_omp/modules"
                elpa_include_omp="$elpa_include_omp/modules"
                ELPA_CFLAGS_OMP="-I'$elpa_include_omp'"
            else
                echo "Cannot find elpa_openmp-${elpa_ver} from paths $INCLUDE_PATHS"
                exit 1
            fi
        fi
        # get the lib paths
        add_lib_from_paths ELPA_LDFLAGS "libelpa.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking ELPA to user paths ===================="
        pkg_install_dir="$with_elpa"
        check_dir "${pkg_install_dir}/include"
        check_dir "${pkg_install_dir}/lib"
        user_include_path="$pkg_install_dir/include"
        elpa_include="$(find_in_paths "elpa-*" user_include_path)"
        if [ "$elpa_include" != "__FALSE__" ] ; then
            echo "ELPA include directory is found to be $elpa_include/modules"
            elpa_include="$elpa_include/modules"
            check_dir "$elpa_include"
            ELPA_CFLAGS="-I'$elpa_include'"
        else
            echo "Cannot find elpa-* from path $user_include_path"
            exit 1
        fi
        if [ "$ENABLE_OMP" = "__TRUE__" ] ; then
            elpa_include_omp="$(find_in_paths "elpa_openmp-*" user_include_path)"
            if [ "$elpa_include_omp" != "__FALSE__" ] ; then
                echo "ELPA include directory threaded version is found to be $elpa_include_omp/modules"
                elpa_include_omp="$elpa_include_omp/modules"
                check_dir "$elpa_include_omp"
                ELPA_CFLAGS_OMP="-I'$elpa_include_omp'"
            else
                echo "Cannot find elpa_openmp-* from path $user_include_path"
                exit 1
            fi
        fi
        ELPA_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        ;;
esac
if [ "$with_elpa" != "__DONTUSE__" ] ; then
    ELPA_LIBS="-lelpa"
    ELPA_LIBS_OMP="-lelpa_openmp"
    cat <<EOF > "${BUILDDIR}/setup_elpa"
prepend_path CPATH "$elpa_include"
prepend_path CPATH "$elpa_include_omp"
EOF
    if [ "$with_elpa" != "__SYSTEM__" ] ; then
        cat <<EOF >> "${BUILDDIR}/setup_elpa"
prepend_path PATH "$pkg_install_dir/bin"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
EOF
    fi
    cat "${BUILDDIR}/setup_elpa" >> $SETUPFILE
    cat <<EOF >> "${BUILDDIR}/setup_elpa"
export ELPA_CFLAGS="${ELPA_CFLAGS}"
export ELPA_LDFLAGS="${ELPA_LDFLAGS}"
export ELPA_LIBS="${ELPA_LIBS}"
export ELPA_CFLAGS_OMP="${ELPA_CFLAGS_OMP}"
export ELPA_LDFLAGS_OMP="${ELPA_LDFLAGS_OMP}"
export ELPA_LIBS_OMP="${ELPA_LIBS_OMP}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(-D__ELPA3|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(IF_OMP(${ELPA_CFLAGS_OMP}|${ELPA_CFLAGS})|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(${ELPA_LDFLAGS}|)"
export CP_LIBS="IF_MPI(IF_OMP(${ELPA_LIBS_OMP}|${ELPA_LIBS})|) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
