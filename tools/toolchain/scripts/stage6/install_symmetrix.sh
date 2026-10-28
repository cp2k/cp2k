#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# symmetrix vendors Kokkos and sphericart as git submodules, so the tarball
# hosted on cp2k.org is a git archive of the pinned commit plus every submodule.
# This builds the host (CPU) libsymmetrix; a Kokkos/GPU build needs CUDA and is
# left as a manual step (see docs/methods/machine_learning/symmetrix.md).
#
# Bump the version, the tarball and the sha256 together with the backend when a
# newer symmetrix is required.
symmetrix_ver="2026.01"
symmetrix_dir="symmetrix-${symmetrix_ver}"
symmetrix_pkg="symmetrix-${symmetrix_ver}.tar.gz"
symmetrix_sha256="56e46d8ebeed4635907c32973a2f1fab2004195d7d8452d120ded074727f3f7f"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_symmetrix" ] && rm "${BUILDDIR}/setup_symmetrix"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "$with_symmetrix" in
  __INSTALL__)
    echo "==================== Installing symmetrix ===================="
    pkg_install_dir="${INSTALLDIR}/${symmetrix_dir}"
    install_lock_file="${pkg_install_dir}/install_successful"
    symmetrix_root="${pkg_install_dir}/${symmetrix_dir}"
    if verify_checksums "${install_lock_file}"; then
      echo "${symmetrix_dir} is already installed, skipping it."
    else
      retrieve_package "${symmetrix_sha256}" "${symmetrix_pkg}"
      echo "Installing from scratch into ${pkg_install_dir}"
      [ -d "${pkg_install_dir}" ] && rm -rf "${pkg_install_dir}"
      mkdir -p "${pkg_install_dir}"
      tar -xzf "${symmetrix_pkg}" -C "${pkg_install_dir}"
      cd "${symmetrix_root}"

      # symmetrix hard-codes "-march=native" as a PRIVATE compile option (its
      # own "TODO: generalize"), which lands after CMAKE_CXX_FLAGS and so wins
      # over anything passed in. Left alone it silently ignores --target-cpu and
      # yields a library that dies with SIGILL on any host older than the build
      # machine. Honour the toolchain's target instead, using -mtune like the
      # other packages do.
      if [ "${TARGET_CPU}" != "native" ]; then
        sed -i "s|-march=native|-mtune=${TARGET_CPU}|" libsymmetrix/CMakeLists.txt
      fi

      # Host (CPU) build; symmetrix requires C++20 and libstdc++'s special math
      # functions (std::assoc_legendre), both provided by the toolchain gcc.
      cmake -S libsymmetrix -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DSYMMETRIX_KOKKOS=OFF \
        > cmake-symmetrix.log 2>&1 || tail_excerpt cmake-symmetrix.log
      cmake --build build -j ${NPROCS:-16} \
        > build-symmetrix.log 2>&1 || tail_excerpt build-symmetrix.log

      cd "${BUILDDIR}"
      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding symmetrix from system paths ===================="
    report_error "symmetrix has no system package; point --with-symmetrix at a prebuilt checkout instead."
    exit 1
    ;;
  __DONTUSE__) ;;
  *)
    echo "==================== Linking symmetrix to user paths ===================="
    symmetrix_root="$with_symmetrix"
    check_dir "${symmetrix_root}/libsymmetrix"
    check_dir "${symmetrix_root}/build"
    ;;
esac

if [ "$with_symmetrix" != "__DONTUSE__" ]; then
  # FindSymmetrix derives the include dirs and static libraries from the
  # checkout + build tree, so the toolchain exports the roots rather than a
  # lib/include prefix. Configure CP2K with -DCP2K_USE_SYMMETRIX=ON; ROOT and
  # BUILD_DIR are picked up from the environment (see FindSymmetrix.cmake). The
  # __SYMMETRIX define is added by the CMake build, not here.
  cat << EOF > "${BUILDDIR}/setup_symmetrix"
export CP2K_SYMMETRIX_ROOT="${symmetrix_root}"
export CP2K_SYMMETRIX_BUILD_DIR="${symmetrix_root}/build"
EOF
  filter_setup "${BUILDDIR}/setup_symmetrix" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_symmetrix"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "symmetrix"

#EOF
