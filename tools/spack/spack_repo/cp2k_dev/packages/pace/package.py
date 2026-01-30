# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Pace(CMakePackage):
    """interatomic Potentials in Atomic Cluster Expansion (PACE)

    This library is required to build the ML-PACE module
    in LAMMPS.

    The library was developed at the Interdisciplinary Centre
    for Advanced Materials Simulation (ICAMS), Ruhr University Bochum

    See `Phys Rev Mat 6 013804 (2022)<https://doi.org/10.1103/PhysRevMaterials.6.013804>`__ and
    `Phys Rev B 99 014104 (2019)<https://doi.org/10.1103/PhysRevB.99.014104>`__ for details.
    """

    maintainers("hjjvandam", "rbberger", "mkrack")

    homepage = (
        "https://www.icams.de/institute/departments-groups/atomistic-modelling-and-simulation/"
    )
    git = "https://github.com/ICAMS/lammps-user-pace.git"

    license("GPL-2.0-or-later", checked_by="hjjvandam")
    version("main", branch="main")
    version("2025.12.4", tag="v.2025.12.4", commit="8d0483447e94d17e01cbc1826a88cae1d8c3aab6")
    version(
        "2023.11.25.2", tag="v.2023.11.25.fix2", commit="e60e850359b918ca93a5e9329548a58d31f4b12b"
    )

    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=False, description="Build shared libraries (otherwise static)")

    depends_on("cmake@3.10:", type="build")
    depends_on("cxx", type="build")

    # Build always with the yaml-cpp and cnpy version included in pace
    depends_on("yaml-cpp")
    depends_on("cnpy")

    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        return args

    def install(self, spec, prefix):
        src = self.stage.source_path
        install_tree(join_path(src, "ML-PACE/ace"), join_path(prefix.include, "ace"))
        install_tree(
            join_path(src, "ML-PACE/ace-evaluator"), join_path(prefix.include, "ace-evaluator")
        )
        install_tree(join_path(src, "yaml-cpp/include"), prefix.include)
        mkdirp(prefix.lib)
        install(join_path(self.build_directory, "libpace.a"), prefix.lib)
