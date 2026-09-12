#!/bin/bash -e

# author: Thomas D. Kuehne, tkuehne@cp2k.org
# shellcheck disable=SC1091

SCRIPT_NAME="${BASH_SOURCE[0]}"
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_NAME}")/.." && pwd -P)"
wannier90_ver="4.0.2"
wannier90_sha256="2d48b371eefa8b58a6c8088c1bdffc13fe3e761111e15c8566e2ee055d8bcdb0"

source "${SCRIPT_DIR}/common_vars.sh"
source "${SCRIPT_DIR}/tool_kit.sh"
source "${SCRIPT_DIR}/signal_trap.sh"
source "${INSTALLDIR}/toolchain.conf"
source "${INSTALLDIR}/toolchain.env"
with_wannier90="${with_wannier90:?Run install_cp2k_toolchain.sh first}"

mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"
if [ -f setup_wannier90 ]; then
  rm setup_wannier90
fi

if [ "${MPI_MODE}" = "no" ]; then
  wannier90_mpi=OFF
  wannier90_pkgconfig=wannier90
else
  wannier90_mpi=ON
  wannier90_pkgconfig=wannier90_mpi
fi

case "${with_wannier90}" in
  __INSTALL__)
    echo "==================== Installing Wannier90 ===================="
    pkg_install_dir="${INSTALLDIR}/wannier90-${wannier90_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"
    if verify_checksums "${install_lock_file}"; then
      echo "Wannier90 is already installed, skipping it."
    else
      archive="wannier90-${wannier90_ver}.tar.gz"
      if ! checksum "${wannier90_sha256}" "${archive}"; then
        download_pkg_from_urlpath "${wannier90_sha256}" "v${wannier90_ver}.tar.gz" \
          "https://github.com/wannier-developers/wannier90/archive/refs/tags" "${archive}"
      fi
      # Examples and reference data are large and are not needed for the library build.
      tar -xzf "${archive}" "wannier90-${wannier90_ver}/CMakeLists.txt" \
        "wannier90-${wannier90_ver}/LICENSE" \
        "wannier90-${wannier90_ver}/cmake" "wannier90-${wannier90_ver}/src"
      cd "wannier90-${wannier90_ver}"
      # toolchain.conf is part of the lock: rebuild if the compiler or MPI mode changed.
      rm -rf build "${pkg_install_dir}"
      mkdir build
      case "${MATH_MODE}" in
        openblas) blas_vendor=OpenBLAS ;;
        mkl) blas_vendor=Intel10_64lp ;;
        acml) blas_vendor=ACML ;;
        *) blas_vendor=All ;;
      esac
      cmake -S . -B build \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBLA_SIZEOF_INTEGER=4 \
        -DBLA_VENDOR="${blas_vendor}" \
        -DWANNIER90_MPI="${wannier90_mpi}" \
        -DWANNIER90_SHARED_LIBS=OFF \
        -DWANNIER90_INSTALL=ON \
        -DWANNIER90_TEST=OFF \
        > build/configure.log 2>&1 || tail_excerpt build/configure.log
      cmake --build build --parallel "$(get_nprocs)" \
        > build/build.log 2>&1 || tail_excerpt build/build.log
      cmake --install build \
        > build/install.log 2>&1 || tail_excerpt build/install.log
      mkdir -p "${pkg_install_dir}/share/licenses/wannier90"
      cp LICENSE "${pkg_install_dir}/share/licenses/wannier90/LICENSE"
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage7/install_wannier90.sh" \
        "${INSTALLDIR}/toolchain.conf"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding Wannier90 from system paths ===================="
    check_pkgconfig "${wannier90_pkgconfig}"
    # The v4.0.2 .pc template uses colons instead of assignments for prefix/libdir.
    # pcfiledir is supplied by pkg-config itself and also handles multiarch libdirs.
    lib_dir="$(dirname "$(pkg-config --variable=pcfiledir "${wannier90_pkgconfig}")")"
    pkg_install_dir="$(pkg-config --variable=prefix "${wannier90_pkgconfig}")"
    if [ -z "${pkg_install_dir}" ]; then
      case "${lib_dir}" in
        */lib/*) pkg_install_dir="${lib_dir%/lib/*}" ;;
        *) pkg_install_dir="$(dirname "${lib_dir}")" ;;
      esac
    fi
    check_dir "${lib_dir}/cmake/Wannier90"
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking Wannier90 to user paths ===================="
    pkg_install_dir="${with_wannier90}"
    ;;
esac

if [ "${with_wannier90}" != "__DONTUSE__" ]; then
  if [ "${with_wannier90}" != "__SYSTEM__" ]; then
    lib_dir="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && lib_dir="${pkg_install_dir}/lib64"
    check_dir "${lib_dir}/cmake/Wannier90"
    cat << EOF > "${BUILDDIR}/setup_wannier90"
prepend_path LD_LIBRARY_PATH "${lib_dir}"
prepend_path LD_RUN_PATH "${lib_dir}"
prepend_path LIBRARY_PATH "${lib_dir}"
prepend_path PKG_CONFIG_PATH "${lib_dir}/pkgconfig"
prepend_path PATH "${pkg_install_dir}/bin"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  cat << EOF >> "${BUILDDIR}/setup_wannier90"
export WANNIER90_ROOT="${pkg_install_dir}"
export WANNIER90_VER="${wannier90_ver}"
EOF
  filter_setup "${BUILDDIR}/setup_wannier90" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_wannier90"
write_toolchain_env "${INSTALLDIR}"
cd "${ROOTDIR}"
report_timing "wannier90"
