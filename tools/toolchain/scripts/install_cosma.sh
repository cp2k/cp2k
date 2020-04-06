#!/bin/bash -e
[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")" && pwd -P)"

cosma_ver="2.0.6"
cosma_sha256="25c73035628f5652a4df99881993c51af5806c9a7ecb27a9b0b9bea78c4f45f4"
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_cosma" ] && rm "${BUILDDIR}/setup_cosma"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_cosma" in
    __INSTALL__)
	require_env OPENBLASROOT
	require_env SCALAPACKROOT
      	
	echo "==================== Installing cosma ===================="
        pkg_install_dir="${INSTALLDIR}/cosma-${cosma_ver}"
        install_lock_file="$pkg_install_dir/install_successful"
        if verify_checksums "${install_lock_file}" ; then
            echo "cosma-${cosma_ver} is already installed, skipping it."
        else
            if [ -f cosma-${cosma_ver}.tar.gz ] ; then
                echo "cosma-${cosma_ver}.tar.gz is found"
	    else
                download_pkg ${DOWNLOADER_FLAGS} ${cosma_sha256} \
                             "https://github.com/eth-cscs/COSMA/releases/download/v${cosma_ver}/cosma.tar.gz" \
                              -o cosma-${cosma_ver}.tar.gz

            fi
            echo "Installing from scratch into ${pkg_install_dir}"
            [ -d cosma-${cosma_ver} ] && rm -rf cosma-${cosma_ver}
            tar -xzf cosma-${cosma_ver}.tar.gz
            mv cosma cosma-${cosma_ver}
            cd cosma-${cosma_ver}
            mkdir build-cpu
            cd build-cpu
            case "$FAST_MATH_MODE" in
                mkl)
                cmake -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
                          -DCOSMA_BLAS=MKL \
                          -DCOSMA_SCALAPACK=MKL \
                          -DCOSMA_WITH_TESTS=NO \
			  -DCOSMA_WITH_APPS=NO \
			  -DCOSMA_WITH_BENCHMARKS=NO .. > cmake.log 2>&1
                    ;;
                *)
                    cmake -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
                          -DCOSMA_BLAS=OPENBLAS \
                          -DCOSMA_SCALAPACK=CUSTOM \
			  -DCOSMA_WITH_TESTS=NO \
			  -DCOSMA_WITH_APPS=NO \
			  -DCMAKE_WITH_BENCHMARKS=NO .. > cmake.log 2>&1
                    ;;
            esac

            make -j $NPROCS > make.log 2>&1
            make -j $NPROCS install > install.log 2>&1
            cd ..
            write_checksums "${install_lock_file}" "${SCRIPT_DIR}/$(basename ${SCRIPT_NAME})"
            if [ "$ENABLE_CUDA" = "__TRUE__" ] ; then
                [ -d build-cuda ] && rm -rf "build-cuda"
                mkdir build-cuda
                cd build-cuda
                case "$FAST_MATH_MODE" in
                 mkl)
                 cmake -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
 			 -DCOSMA_BLAS=CUDA \
                         -DCOSMA_SCALAPACK=MKL \
			 -DCOSMA_WITH_TESTS=NO \
			 -DCOSMA_WITH_BENCHMARKS=NO \
			 -DCOSMA_WITH_APPS=NO .. > cmake.log 2>&1
                     ;;
                 *)
                     cmake -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
                           -DCOSMA_BLAS=CUDA \
                           -DCOSMA_SCALAPACK=CUSTOM \
                           -DCOSMA_WITH_TESTS=NO \
			   -DCOSMA_WITH_BENCHMARKS=NO \
			   -DCOSMA_WITH_APPS=NO .. > cmake.log 2>&1
                     ;;
             esac
                make -j $NPROCS > make.log 2>&1
                install -d ${pkg_install_dir}/lib/cuda
                [ -f libs/grid2grid/src/grid2grid/*.a ] && install -m 644 libs/grid2grid/src/grid2grid/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
                [ -f libs/options/*.so ] && install -m 644 libs/options/*.so ${pkg_install_dir}/lib/cuda >> install.log 2>&1
                install -m 644 src/cosma/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
		[ -f libs/Tiled-MM/src/Tiled-MM/*.a ] && install -m 644 libs/Tiled-MM/src/Tiled-MM/*.a ${pkg_install_dir}/lib/cuda >> install.log 2>&1
            fi
        fi
        COSMA_ROOT="${pkg_install_dir}"
        COSMA_CFLAGS="-I'${pkg_install_dir}/include'"

        # check if cosma is compiled with 64bits and set up COSMA_LIBDIR accordingly
        COSMA_LIBDIR="${pkg_install_dir}/lib"

        [ -d  ${pkg_install_dir}/lib64 ] && COSMA_LIBDIR="${pkg_install_dir}/lib64"

        COSMA_LDFLAGS="-L'${COSMA_LIBDIR}' -Wl,-rpath='${COSMA_LIBDIR}'"
        COSMA_CUDA_LDFLAGS="-L'${COSMA_LIBDIR}/cuda' -Wl,-rpath='${COSMA_LIBDIR}/cuda'"
        ;;
    __SYSTEM__)
        echo "==================== Finding cosma from system paths ===================="
        check_command pkg-config --modversion cosma
        add_include_from_paths COSMA_CFLAGS "cosma.h" $INCLUDE_PATHS
        add_lib_from_paths COSMA_LDFLAGS "libcosma.*" $LIB_PATHS
        ;;
    __DONTUSE__)
        ;;
    *)
        echo "==================== Linking cosma to user paths ===================="
        pkg_install_dir="$with_cosma"
        check_dir "$pkg_install_dir/lib"
        check_dir "$pkg_install_dir/lib64"
        check_dir "$pkg_install_dir/include"


        # check if cosma is compiled with 64bits and set up COSMA_LIBDIR accordingly
        COSMA_LIBDIR="${pkg_install_dir}/lib"

        [ -d  ${pkg_install_dir}/lib64 ] && COSMA_LIBDIR="${pkg_install_dir}/lib64"

        COSMA_CFLAGS="-I'${pkg_install_dir}/include'"
        COSMA_LDFLAGS="-L'${COSMA_LIBDIR}' -Wl,-rpath='${COSMA_LIBDIR}'"
        ;;
esac
    if [ "$with_cosma" != "__SYSTEM__" ] ; then
        cat << EOF > "${BUILDDIR}/setup_cosma"
prepend_path LD_LIBRARY_PATH "${COSMA_LIBDIR}"
prepend_path LD_RUN_PATH "${COSMA_LIBDIR}"
prepend_path LIBRARY_PATH "${COSMA_LIBDIR}"
prepend_path CPATH "$pkg_install_dir/include"
export COSMA_INCLUDE_DIR="$pkg_install_dir/include"
export COSMA_LIBS="-lcosma_pxgemm -lcosma -lgrid2grid IF_CUDA(-lTiled-MM|)"
export COSMA_ROOT="${pkg_install_dir}"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${COSMA_LIBDIR}/pkgconfig"
EOF
    fi
    cat << EOF >> "${BUILDDIR}/setup_cosma"
export COSMA_GPU_TILE_M=128
export COSMA_GPU_TILE_N=128
export COSMA_GPU_TILE_K=128
export COSMA_CFLAGS="${COSMA_CFLAGS}"
export COSMA_LDFLAGS="${COSMA_LDFLAGS}"
export COSMA_CUDA_LDFLAGS="${COSMA_CUDA_LDFLAGS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_MPI(IF_OMP(-D__COSMA|)|)"
export CP_CFLAGS="\${CP_CFLAGS} ${COSMA_CFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} IF_CUDA(${COSMA_CUDA_LDFLAGS}|${COSMA_LDFLAGS})"
export COSMA_LIBRARY="-lcosma_pxgemm -lcosma -lgrid2grid IF_CUDA(-lTiled-MM|)"
export COSMA_ROOT="$pkg_install_dir"
export COSMA_INCLUDE_DIR="$pkg_install_dir/include"
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:${COSMA_LIBDIR}/pkgconfig"
export COSMA_VERSION=${cosma_ver}
export CP_LIBS="IF_MPI(IF_OMP(${COSMA_LIBS}|)|) \${CP_LIBS}"
EOF
    cat "${BUILDDIR}/setup_cosma" >> $SETUPFILE
fi

# update toolchain environment
load "${BUILDDIR}/setup_cosma"
export -p > "${INSTALLDIR}/toolchain.env"

cd "${ROOTDIR}"
report_timing "cosma"
