#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

sirius_ver=${sirius_ver:-5.8.3}
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh

with_sirius=${1:-__INSTALL__}

if [ $MPI_MODE = no ] ; then
    if [ "$ENABLE_OMP" = "__FALSE__"] ; then
        report_warning $LINENO "MPI is disabled, skipping sirius installation"
cat <<EOF > "${BUILDDIR}/setup_sirius"
    with_sirius="__FALSE__"
EOF
      exit 0
    fi
fi
[ -f "${BUILDDIR}/setup_sirius" ] && rm "${BUILDDIR}/setup_sirius"

SIRIUS_CFLAGS=''
SIRIUS_LDFLAGS=''
SIRIUS_LIBS=''

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
case "$with_sirius" in
    __INSTALL__)
        echo "==================== Installing SIRIUS ========================================="
        echo "                                                                                "
        echo " note that this package depends on gsl, libspg, elpa, scalapack, json-fortran   "
        echo " hdf5 and libxc. this package can not work without these libraries dependencies "
        echo " activated.                                                                     "
        echo "                                                                                "
        echo "================================================================================"


        require_env FFTW_LDFLAGS
        require_env FFTW_LIBS
        require_env FFTW_CFLAGS
        require_env ELPA_LDFLAGS
        require_env ELPA_LIBS
        require_env ELPHA_CFLAGS
        require_env GSL_LDFLAGS
        require_env GSL_CFLAGS
        require_env GSL_LIBS
        require_env MATH_LIBS
        require_env MPI_LDFLAGS
        require_env MPI_LIBS
        require_env SCALAPACK_LDFLAGS
        require_env SCALAPACK_CFLAGS
        require_env SCALAPACK_LIBS
        require_env LIBXC_LIBS
        require_env LIBXC_CFLAGS
        require_env LIBXC_LDFLAGS
        require_env SPGLIB_LIBS
        require_env SPGLIB_CFLAGS
        require_env SPGLIB_LDFLAGS
        require_env HDF5_LIBS
        require_env HDF5_CFLAGS
        require_env HDF5_LDFLAGS

        pkg_install_dir="${INSTALLDIR}/sirius-${sirius_ver}"
        install_lock_file="${pkg_install_dir}/install_successful"
        if [ -f "${install_lock_file}" ] ; then
            echo "sirius_dist-${sirius_ver} is already installed, skipping it."
        else
            if [ -f SIRIUS-${sirius_ver}.tar.gz ] ; then
                echo "sirius_v${sirius_ver}.tar.gz is found"
            else
                download_pkg -o SIRIUS-${sirius_ver}.tar.gz \
                             https://github.com/electronic-structure/SIRIUS/archive/v${sirius_ver}.tar.gz
            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d sirius-${sirius_ver} ] && rm -rf sirius-${sirius_ver}
            tar -xzf SIRIUS-${sirius_ver}.tar.gz
            cd SIRIUS-${sirius_ver}
            # rm -fr build
            # mkdir build
            # cd build
            # if [ "$ENABLE_CUDA" = __TRUE__ ] ; then
            #     cmake -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} ..
            # else
            #     cmake -DCMAKE_INSTALL_PREFIX=${pkg_install_dir} -DUSE_CUDA=off ..
            # fi
            # make
            # make install
            # cd ../..
            #touch "${install_lock_file}"
            cat <<EOF > make.inc
CXX = \$(MPICXX)
BASIC_CXX_OPT = -O3 -DNDEBUG -mtune=native -ftree-loop-vectorize ${MATH_CFLAGS}
CXX_OPT = \$(BASIC_CXX_OPT) -fopenmp -std=c++11 -D__SCALAPACK -D__ELPA
EOF
             if [ "$ENABLE_CUDA" = __TRUE__ ] ; then
                 cat <<EOF >> make.inc
CXX_OPT := \$(CXX_OPT) -D__GPU -I${CUDA_SDK}/include -I${CUDA_ROOT}/include
NVCC=nvcc -O3 -arch=sm_60
LIBS := ${CUDA_LIBS}
CXX_OPT := \$(CXX_OPT) ${CUDA_CFLAGS}
EOF
            fi
            cat <<EOF >> make.inc
CXX_OPT := \$(CXX_OPT) -I${PWD}/src
CXX_OPT := \$(CXX_OPT) -I${PWD}/src/SDDK
CXX_OPT := \$(CXX_OPT) ${ELPA_CFLAGS}
CXX_OPT := \$(CXX_OPT) ${GSL_CFLAGS}
CXX_OPT := \$(CXX_OPT) ${SPGLIB_CFLAGS}
CXX_OPT := \$(CXX_OPT) ${HDF5_CFLAGS}
MPI_FC = ${MPIFC}
MPI_FC_OPT = -g -O2 -fopenmp -cpp

LIBS := ${LIBXC_LDFLAGS} ${FFTW_LDFLAGS} ${ELPA_LDFLAGS} ${SCALAPACK_LDFLAGS} ${GSL_LDFLAGS} ${SPG_LDFLAGS}
LIBS := \$(LIBS) ${SPGLIB_LIBS} ${FFTW_LIBS} ${GSL_LIBS} ${SCALAPACK_LIBS} ${ELPA_LIBS}
EOF

            # this a hack needs to be fixed permanently
            sed -i -e "s/: log/: /g" src/Makefile
cat > src/version.hpp <<EOF
#ifndef __VERSION_H__
#define __VERSION_H__
const char* const git_hash = "6ec392682d7ffda6090ee3124d6841adc6a5fd0d";
const char* const git_branchname = "heads/develop";
const char* const build_date = "Thu, 23 Aug 2018 13:58:38";
#endif
EOF

             # a fix to indicate that spglib has changed the directory where the header file is located
             sed -i -e "s/spglib\///g" src/Unit_cell/unit_cell_symmetry.hpp

             # small patch to fix an error in sirius. It is already fixed in the develop branch. just waiting for a new release
             sed -i -e "s/mpi_fin = call_mpi_fin/mpi_fin = *call_mpi_fin/g" src/sirius_api.cpp
             sed -i -e "s/device_reset = call_device_reset/device_reset = *call_device_reset/g" src/sirius_api.cpp
             sed -i -e "s/fftw_fin = call_fftw_fin/fftw_fin = *call_fftw_fin/g" src/sirius_api.cpp
             make -C src
             install -v -d ${pkg_install_dir}/include
             install -v -d ${pkg_install_dir}/lib
             cp -Rv src/* ${pkg_install_dir}/include
             rm -f ${pkg_install_dir}/include/*.f90
             #rm -f ${pkg_install_dir}/include/*.mod
             rm -f ${pkg_install_dir}/include/*.o
             install -v -m 644 src/*.a ${pkg_install_dir}/lib
             install -v -m 644 src/*.mod ${pkg_install_dir}/lib
        fi
        SIRIUS_CFLAGS="-I'${pkg_install_dir}/include'"
        SIRIUS_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
        touch "${install_lock_file}"
        ;;
    __SYSTEM__)
        echo "==================== Finding SIRIUS_DIST from system paths ===================="
        echo "Warning very unlickely you have this package installed by default              "
        echo "                                                                               "
        echo " note that this package depends on gsl, libspg, elpa, scalapack, json-fortran  "
        echo " hdf5 and libxc. this package can not work without these libraries dependencies     "
        echo " activated.                                                                    "
        echo "                                                                               "
        echo "==============================================================================="

        require_env FFTW_LDFLAGS
        require_env FFTW_LIBS
        require_env FFTW_CFLAGS
        require_env ELPA_LDFLAGS
        require_env ELPA_LIBS
        require_env ELPA_CFLAGS
        require_env GSL_LDFLAGS
        require_env GSL_CFLAGS
        require_env GSL_LIBS
        require_env MATH_LIBS
        require_env MPI_LDFLAGS
        require_env MPI_LIBS
        require_env SCALAPACK_LDFLAGS
        require_env SCALAPACK_CFLAGS
        require_env SCALAPACK_LIBS
        require_env LIBXC_LIBS
        require_env LIBXC_CFLAGS
        require_env LIBXC_LDFLAGS
        require_env SPGLIB_LIBS
        require_env SPGLIB_CFLAGS
        require_env SPGLIB_LDFLAGS
        require_env HDF5_LIBS
        require_env HDF5_CFLAGS
        require_env HDF5_LDFLAGS

        check_lib -lsirius "sirius"
        add_include_from_paths SIRIUS_CFLAGS "sirius*" $INCLUDE_PATHS
        add_lib_from_paths SIRIUS_LDFLAGS "libsirius.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking SIRIUS_Dist to user paths ===================="
        pkg_install_dir="$with_sirius"
        check_dir "${pkg_install_dir}/lib"
        check_dir "${pkg_install_dir}/include"
        ;;
esac
if [ "$with_sirius" != "__DONTUSE__" ] ; then
    SIRIUS_LIBS="-lsirius_f -lsirius"
    if [ "$with_sirius" != "__SYSTEM__" ] ; then
        cat <<EOF > "${BUILDDIR}/setup_sirius"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path CPATH "$pkg_install_dir/include"
EOF
        cat "${BUILDDIR}/setup_sirius" >> $SETUPFILE
    fi
    cat <<EOF >> "${BUILDDIR}/setup_sirius"
export SIRIUS_CFLAGS="-I'${pkg_install_dir}/include'"
export SIRIUS_FFLAGS="-I'${pkg_install_dir}/include'"
export SIRIUS_LDFLAGS="${SIRIUS_LDFLAGS}"
export SIRIUS_LIBS="${SIRIUS_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(IF_OMP(-D__SIRIUS|)|)"
export CP_CFLAGS="\${CP_CFLAGS} IF_MPI(IF_OMP(-I'${pkg_install_dir}/include'|)|)"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_MPI(IF_OMP(${SIRIUS_LDFLAGS}|)|)"
export CP_LIBS="IF_MPI(IF_OMP(${SIRIUS_LIBS}|)|) \${CP_LIBS}"
EOF
fi
cd "${ROOTDIR}"
