#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=SC1003,SC1035,SC1083,SC1090
# shellcheck disable=SC2001,SC2002,SC2005,SC2016,SC2091,SC2034,SC2046,SC2086,SC2089,SC2090
# shellcheck disable=SC2124,SC2129,SC2144,SC2153,SC2154,SC2155,SC2163,SC2164,SC2166
# shellcheck disable=SC2235,SC2237

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libsmm" ] && rm "${BUILDDIR}/setup_libsmm"

LIBSMM_CFLAGS=''
LIBSMM_LDFLAGS=''
LIBSMM_LIBS=''
! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_libsmm" in
  __INSTALL__)
    echo "==================== Installing libsmm ===================="
    pkg_install_dir="${INSTALLDIR}/libsmm"
    install_lock_file="$pkg_install_dir/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "libsmm is already installed, skipping it."
    else
      # Here we attempt to determine which precompiled libsmm binary
      # to download, and do that if such binary exists on CP2K web
      # repository.  The binary is determined by the arch and
      # libcore values obtained via OpenBLAS prebuild.
      echo "Searching for an optimised libsmm binary from CP2K website"
      case ${OPENBLAS_LIBCORE} in
        haswell)
          libsmm="libsmm_dnn_haswell-2015-11-10.a"
          libsmm_sha256="a1cf9eb1bfb1bd3467024e47173a6e85881c3908961b7bb29f3348af9837018b"
          echo "An optimized libsmm $libsmm is available"
          ;;
        ivybridge)
          libsmm="libsmm_dnn_ivybridge-2015-07-02.a"
          libsmm_sha256="ef74fb7339979545f9583d9ecab52c640c4f98f9dd49f98d2b4580304d5fcf60"
          echo "An optimized libsmm $libsmm is available"
          ;;
        nehalem)
          libsmm="libsmm_dnn_nehalem-2015-07-02.a"
          libsmm_sha256="cc7e8c6623055fc6bc032dfda2d08b2201a8d86577ab72c3f66bee9b86cbebe9"
          echo "An optimized libsmm $libsmm is available"
          ;;
        sandybridge)
          libsmm="libsmm_dnn_sandybridge-2015-11-10.a"
          libsmm_sha256="56ffdafa715554ec87f20ff1e0150450209a7635b47b8a5b81970e88ec67687c"
          echo "An optimized libsmm $libsmm is available"
          ;;
        *)
          echo "No optimised binary found ..."
          echo "Searching for a generic libsmm binary from CP2K website"
          if [ "${OPENBLAS_ARCH}" = "x86_64" ]; then
            libsmm="libsmm_dnn_x86_64-latest.a"
            libsmm_sha256="dd58aee2bc5505e23b0761835bf2b9a90e5f050c6708ef68c5028373970673f8"
            echo "A generic libsmm $libsmm is available."
            echo "Consider building and contributing to CP2K an optimized"
            echo "libsmm for your $OPENBLAS_ARCH $OPENBLAS_LIBCORE using"
            echo "the toolkit in tools/build_libsmm provided in cp2k package"
          fi
          ;;
      esac
      # we know what to get, proceed with install
      if [ "x$libsmm" != "x" ]; then
        if [ -f $libsmm ]; then
          echo "$libsmm has already been downloaded."
        else
          download_pkg ${DOWNLOADER_FLAGS} $libsmm_sha256 https://www.cp2k.org/static/downloads/libsmm/$libsmm
        fi
        # install manually
        ! [ -d "${pkg_install_dir}/lib" ] && mkdir -p "${pkg_install_dir}/lib"
        cp $libsmm "${pkg_install_dir}/lib"
        ln -sf "${pkg_install_dir}/lib/$libsmm" "${pkg_install_dir}/lib/libsmm_dnn.a"
      else
        echo "No libsmm is available"
        echo "Consider building an optimized libsmm on your system yourself"
        echo "using the toolkid in tools/build_libsmm provided in cp2k package"
        cat << EOF > "${BUILDDIR}/setup_libsmm"
with_libsmm="__DONTUSE__"
EOF
        exit 0
      fi
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage4/$(basename ${SCRIPT_NAME})"
    fi
    LIBSMM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
  __SYSTEM__)
    echo "==================== Finding Libsmm from system paths ===================="
    check_lib -lsmm_dnn "libsmm"
    add_lib_from_paths LIBSMM_LDFLAGS "libsmm_dnn.*" $LIB_PATHS
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking Libsmm to user paths ===================="
    pkg_install_dir="$with_libsmm"
    check_dir "${pkg_install_dir}/lib"
    LIBSMM_LDFLAGS="-L'${pkg_install_dir}/lib' -Wl,-rpath='${pkg_install_dir}/lib'"
    ;;
esac
if [ "$with_libsmm" != "__DONTUSE__" ]; then
  LIBSMM_LIBS="-lsmm_dnn"
  if [ "$with_libsmm" != "__SYSTEM__" ]; then
    cat << EOF > "${BUILDDIR}/setup_libsmm"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
EOF
    cat "${BUILDDIR}/setup_libsmm" >> $SETUPFILE
  fi
  cat << EOF >> "${BUILDDIR}/setup_libsmm"
export LIBSMM_LDFLAGS="${LIBSMM_LDFLAGS}"
export LIBSMM_LIBS="${LIBSMM_LIBS}"
export CP_DFLAGS="\${CP_DFLAGS} IF_VALGRIND(|-D__HAS_smm_dnn)"
export CP_LDFLAGS="\${CP_LDFLAGS} ${LIBSMM_LDFLAGS}"
export CP_LIBS="IF_VALGRIND(|${LIBSMM_LIBS}) \${CP_LIBS}"
EOF
fi

load "${BUILDDIR}/setup_libsmm"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libsmm"
