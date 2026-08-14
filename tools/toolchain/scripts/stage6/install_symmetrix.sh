#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# symmetrix has no tagged release, no download mirror and vendors Kokkos and
# sphericart as git submodules, so it is fetched by cloning a pinned commit with
# --recursive (unlike the tarball-based packages in this stage). This builds the
# host (CPU) libsymmetrix; a Kokkos/GPU build needs CUDA and is left as a manual
# step (see docs/methods/machine_learning/symmetrix.md).
symmetrix_repo="https://github.com/wcwitt/symmetrix"
# Pinned symmetrix commit validated against the CP2K symmetrix backend. Bump this
# together with the backend when a newer symmetrix is required.
symmetrix_commit="739660acb523bc5f30cd57cee6611161a33e2605"

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
    pkg_install_dir="${INSTALLDIR}/symmetrix-${symmetrix_commit}"
    install_lock_file="${pkg_install_dir}/install_successful"
    symmetrix_root="${pkg_install_dir}/symmetrix"
    if verify_checksums "${install_lock_file}"; then
      echo "symmetrix-${symmetrix_commit} is already installed, skipping it."
    else
      echo "Cloning ${symmetrix_repo} @ ${symmetrix_commit} into ${symmetrix_root}"
      [ -d "${pkg_install_dir}" ] && rm -rf "${pkg_install_dir}"
      mkdir -p "${pkg_install_dir}"
      git clone --recursive "${symmetrix_repo}" "${symmetrix_root}" \
        > clone-symmetrix.log 2>&1 || tail_excerpt clone-symmetrix.log
      cd "${symmetrix_root}"
      git checkout "${symmetrix_commit}" \
        > checkout-symmetrix.log 2>&1 || tail_excerpt checkout-symmetrix.log
      git submodule update --init --recursive \
        > submodule-symmetrix.log 2>&1 || tail_excerpt submodule-symmetrix.log

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
