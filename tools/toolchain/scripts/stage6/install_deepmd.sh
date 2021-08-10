#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

deepmd_ver="2.0.0.b1"
deepmd_sha256="4860d72ac8ecd485a501a51a8b3e85f1a62ecf477b164ad636b5824a63c83750"

case "$DEEPMD_MODE" in
  cpu)
    tfcc_ver="2.3.0-cpu_cudaNone_0"
    tfcc_sha256="b1b6e5c68a58e516178e6e855a70dd020021d1cf22a94ad291bb8e7541ec58db"
    ;;
  cuda)
    tfcc_ver="2.3.0-gpu_cuda10.1_0"
    tfcc_sha256="085ba903c76adbbcb57dcdb474e5da5fd36324ca64ce673acacf122d5a2f2b3d"
    ;;
esac

cub_link="https://github.com/NVlabs/cub/archive/c3cceac115c072fb63df1836ff46d8c60d9eb304.zip"
cub_sha256="8894c68d7549681591c34078dcd40cf34459b6c7d33407f07e2145d2adb683ee"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_deepmd" ] && rm "${BUILDDIR}/setup_deepmd"

DEEPMD_LDFLAGS=''
DEEPMD_LIBS=''
DEEPMD_CXXFLAGS='-std=gnu++11 '

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_tfcc" in
  __DONTUSE__)
    tensorflow_root="$with_deepmd"
    check_dir "${tensorflow_root}/include"
    check_dir "${tensorflow_root}/lib"
    ;;
  __INSTALL__)
    echo "==================== Installing TensorFlow C++ Interface ===================="
    pkg_install_dir="${INSTALLDIR}/deepmd-kit-${deepmd_ver}"
    export DP_VARIANT=${DEEPMD_MODE}
    if [ -f libtensorflow_cc-${tfcc_ver}.tar.bz2 ]; then
      echo "libtensorflow_cc-${tfcc_ver}.tar.bz2 is found"
    else
      download_pkg ${DOWNLOADER_FLAGS} ${deepmd_sha256} \
        https://anaconda.org/deepmodeling/libtensorflow_cc/2.3.0/download/linux-64/libtensorflow_cc-${tfcc_ver}.tar.bz2
    fi
    [ -d libtensorflow_cc-${tfcc_ver} ] && rm -rf libtensorflow_cc-${tfcc_ver}
    echo "Installing from scratch into ${pkg_install_dir}"
    tensorflow_root=${pkg_install_dir}/dp
    mkdir -p ${tensorflow_root}
    tar -xjf libtensorflow_cc-${tfcc_ver}.tar.bz2 -C ${tensorflow_root}
    cd ..
    ;;
  __SYSTEM__)
    check_lib -ltensorflow_cc "DEEPMD"
    check_lib -ltensorflow_framework "DEEPMD"
    add_lib_from_paths DEEPMD_LDFLAGS "libtensorflow*" $LIB_PATHS
    add_include_from_paths -p DEEPMD_CFLAGS "tensorflow" $INCLUDE_PATHS
    add_include_from_paths -p DEEPMD_CXXFLAGS "tensorflow" $INCLUDE_PATHS
    ;;
  *)
    tensorflow_root="$with_tfcc"
    check_dir "${tensorflow_root}/include"
    check_dir "${tensorflow_root}/lib"
    ;;
esac

cd "${BUILDDIR}"
case "$with_deepmd" in
  __INSTALL__)
    echo "==================== Installing DeePMD ===================="
    pkg_install_dir="${INSTALLDIR}/deepmd-kit-${deepmd_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    deepmd_root="${pkg_install_dir}/dp"
    if verify_checksums "${install_lock_file}"; then
      echo "deepmd-kit-${deepmd_ver} is already installed, skipping it."
    else
      if [ -f deepmd-kit-v${deepmd_ver}.tar.gz ]; then
        echo "deepmd-kit-v${deepmd_ver}.tar.gz is found"
      else
        download_pkg ${DOWNLOADER_FLAGS} -o deepmd-kit-v${deepmd_ver}.tar.gz ${deepmd_sha256} \
          https://github.com/deepmodeling/deepmd-kit/archive/refs/tags/v${deepmd_ver}.tar.gz
      fi
      [ -d deepmd-kit-${deepmd_ver} ] && rm -rf deepmd-kit-${deepmd_ver}
      echo "Installing from scratch into ${pkg_install_dir}"
      tar -xzf deepmd-kit-v${deepmd_ver}.tar.gz
      if [ "${DEEPMD_MODE}" == "cuda" ]; then
        if [ -f cub.zip ]; then
          echo "cub.zip is found"
        else
          download_pkg ${DOWNLOADER_FLAGS} -o cub.zip ${cub_sha256} ${cub_link}
        fi
        rm -r ${BUILDDIR}/deepmd-kit-${deepmd_ver}/source/lib/src/cuda/cub
        unzip cub.zip -d ${BUILDDIR}/deepmd-kit-${deepmd_ver}/source/lib/src/cuda
        mv ${BUILDDIR}/deepmd-kit-${deepmd_ver}/source/lib/src/cuda/cub-c3cceac115c072fb63df1836ff46d8c60d9eb304 \
          ${BUILDDIR}/deepmd-kit-${deepmd_ver}/source/lib/src/cuda/cub
      fi
      mkdir -p ${deepmd_root}
      mkdir -p ${BUILDDIR}/deepmd-kit-${deepmd_ver}/source/build
      cd ${BUILDDIR}/deepmd-kit-${deepmd_ver}/source/build
      cmake -DTENSORFLOW_ROOT=${deepmd_root} -DCMAKE_INSTALL_PREFIX=${tensorflow_root} ..
      NPROC=$(nproc --all)
      make -j${NPROC}
      make install
      cd ../../..
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    DEEPMD_DFLAGS="-D__DEEPMD -DHIGH_PREC"
    DEEPMD_CFLAGS="-I'${deepmd_root}/include/deepmd/' -I'${tensorflow_root}/include'"
    DEEPMD_CXXFLAGS="-std=gnu++11 -I'${deepmd_root}/include/deepmd/' -I'${tensorflow_root}/include'"
    DEEPMD_LDFLAGS="-L'${deepmd_root}/lib' -L'${tensorflow_root}/lib' -Wl,--no-as-needed -Wl,-rpath='${deepmd_root}/lib' -Wl,-rpath='${tensorflow_root}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding DeePMD from system paths ===================="
    check_lib -ldeepmd "DEEPMD"
    add_lib_from_paths DEEPMD_LDFLAGS "libdeepmd*" $LIB_PATHS
    add_include_from_paths DEEPMD_CFLAGS "deepmd" $INCLUDE_PATHS
    add_include_from_paths DEEPMD_CXXFLAGS "deepmd" $INCLUDE_PATHS
    check_lib -ltensorflow_cc "DEEPMD"
    check_lib -ltensorflow_framework "DEEPMD"
    add_lib_from_paths DEEPMD_LDFLAGS "libtensorflow*" $LIB_PATHS
    DEEPMD_DFLAGS="-D__DEEPMD -DHIGH_PREC"
    add_include_from_paths -p DEEPMD_CFLAGS "tensorflow" $INCLUDE_PATHS
    add_include_from_paths -p DEEPMD_CXXFLAGS "tensorflow" $INCLUDE_PATHS
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking DEEPMD to user paths ===================="
    deepmd_root="$with_deepmd"
    check_dir "${deepmd_root}/include/deepmd"
    check_dir "${deepmd_root}/lib"
    DEEPMD_DFLAGS="-D__DEEPMD -DHIGH_PREC"
    DEEPMD_CFLAGS="-I'${deepmd_root}/include/deepmd/' -I'${tensorflow_root}/include'"
    DEEPMD_CXXFLAGS="-std=gnu++11 -I'${deepmd_root}/include/deepmd/' -I'${tensorflow_root}/include'"
    DEEPMD_LDFLAGS="-L'${deepmd_root}/lib' -L'${tensorflow_root}/lib' -Wl,--no-as-needed -Wl,-rpath='${deepmd_root}/lib' -Wl,-rpath='${tensorflow_root}/lib'"
    ;;
esac

if [ "$with_deepmd" != "__DONTUSE__" ]; then
  if [ "$DEEPMD_MODE" == "cpu" ]; then
    DEEPMD_LIBS='-ldeepmd_op -ldeepmd -ldeepmd_cc -ltensorflow_cc -ltensorflow_framework -lstdc++'
  elif [ "$DEEPMD_MODE" == "cuda" ]; then
    DEEPMD_LIBS='-ldeepmd_op -ldeepmd -ldeepmd_cc -ldeepmd_op_cuda -ltensorflow_cc -ltensorflow_framework -lstdc++'
  fi
  if [ "$with_deepmd" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_deepmd"
prepend_path LD_LIBRARY_PATH "$pkg_install_dir/lib"
prepend_path LD_RUN_PATH "$pkg_install_dir/lib"
prepend_path LIBRARY_PATH "$pkg_install_dir/lib"
EOF
    cat "${BUILDDIR}/setup_deepmd" >> $SETUPFILE
  fi

  cat << EOF >> "${BUILDDIR}/setup_deepmd"
export DEEPMD_DFLAGS="${DEEPMD_DFLAGS}"
export DEEPMD_CFLAGS="${DEEPMD_CFLAGS}"
export DEEPMD_CXXFLAGS="${DEEPMD_CXXFLAGS}"
export DEEPMD_LDFLAGS="${DEEPMD_LDFLAGS}"
export DEEPMD_LIBS="${DEEPMD_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} ${DEEPMD_DFLAGS}"
export CP_CFLAGS="\${CP_CFLAGS} ${DEEPMD_CFLAGS}"
export CP_CXXFLAGS="\${CP_CXXFLAGS} ${DEEPMD_CXXFLAGS}"
export CP_LDFLAGS="\${CP_LDFLAGS} ${DEEPMD_LDFLAGS}"
export CP_LIBS="\${CP_LIBS} ${DEEPMD_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_deepmd"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "deepmd"
