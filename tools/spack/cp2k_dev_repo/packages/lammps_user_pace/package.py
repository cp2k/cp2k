# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.cmake import CMakePackage


class LammpsUserPace(CMakePackage):
    """Atomic Cluster Expansion (ACE) for carbon. The ACE was parameterized based
    on an exhaustive set of important carbon structures that were computed using
    accurate density functional theory (DFT) calculations."""

    homepage = "https://github.com/ICAMS/lammps-user-pace"
    url = "https://github.com/ICAMS/lammps-user-pace/archive/v.2023.11.25.fix2.tar.gz"
    git = "https://github.com/ICAMS/lammps-user-pace.git"

    license("GPL-2.0-only OR GPL-3.0-only", checked_by="mkrack")

    maintainers("mkrack")

    version(
        "2023.11.25.fix2",
        sha256="e0885351a8a730f5576dace2374fa470523a4526383c6a64af571e1344a40686",
    )

    variant("cxx_standard", default="11", description="Required CXX standard")
    variant("pic", default=True, description="Compile position independent code (PIC)")

    depends_on("cxx", type="build")

    def cmake_args(self):
        args = [
            self.define_from_variant("CXX_STANDARD", "cxx_standard"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        return args

    def install(self, spec, prefix):
        mkdirp(join_path(prefix.include, "ace"))
        install("ML-PACE/ace/*.h", join_path(prefix.include, "ace"))
        mkdirp(join_path(prefix.include, "ace-evaluator"))
        install("ML-PACE/ace-evaluator/*.h", join_path(prefix.include, "ace-evaluator"))
        copy_tree("yaml-cpp/include", prefix.include)
        with working_dir(prefix.lib, create=True):
            for lib in ["libpace.a", "libcnpy.a", "build-yaml-cpp/libyaml-cpp-pace.a"]:
                install(join_path(self.build_directory, lib), prefix.lib)
