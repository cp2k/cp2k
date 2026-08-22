# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Symmetrix(CMakePackage):
    """Torch-free C++ (optionally Kokkos/GPU) evaluator for MACE interatomic
    potentials.

    symmetrix links against a small static library and lets CP2K run MACE models
    as a FIST many-body potential without depending on LibTorch. It vendors
    Kokkos and sphericart as git submodules and has no install()/package-config
    rules of its own, so this package builds libsymmetrix in place and stages the
    checkout plus its build tree into the prefix. CP2K's FindSymmetrix module
    then derives the include directories and static libraries from that tree
    (point CP2K at it with -DCP2K_SYMMETRIX_ROOT=<prefix>).
    """

    homepage = "https://github.com/wcwitt/symmetrix"
    git = "https://github.com/wcwitt/symmetrix.git"

    maintainers("wcwitt")

    license("MIT")

    # Pinned commit validated against the CP2K symmetrix backend; submodules
    # (Kokkos, sphericart) are fetched because symmetrix has no release tarball.
    # Update the version string and commit together when bumping symmetrix.
    version("2026.01", commit="739660acb523bc5f30cd57cee6611161a33e2605", submodules=True)

    variant(
        "kokkos",
        default=False,
        description="Build the Kokkos (GPU) libsymmetrix evaluator (needs CUDA)",
    )

    # symmetrix requires C++20 and, on the host, a BLAS implementation. CMake's
    # project() enables C as well, so both language dependencies are needed.
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("cmake@3.18:", type="build")
    depends_on("blas")
    depends_on("cuda@12.2:", when="+kokkos")

    # The top-level CMakeLists lives in the libsymmetrix subdirectory.
    root_cmakelists_dir = "libsymmetrix"

    def cmake_args(self):
        return [
            self.define_from_variant("SYMMETRIX_KOKKOS", "kokkos"),
            self.define("CMAKE_CXX_STANDARD", "20"),
        ]

    def install(self, spec, prefix):
        # libsymmetrix has no install target: it is consumed from the checkout +
        # build tree. Preserve that layout in the prefix so FindSymmetrix works
        # with CP2K_SYMMETRIX_ROOT=<prefix> (BUILD_DIR defaults to <prefix>/build).
        install_tree(self.stage.source_path, prefix)
        install_tree(self.build_directory, join_path(prefix, "build"))
