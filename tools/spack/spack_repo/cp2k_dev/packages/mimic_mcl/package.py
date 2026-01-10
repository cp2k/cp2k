# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class MimicMcl(CMakePackage):
    """MiMiC is a high-performance framework for multiscale molecular dynamics simulations.
    MCL, or MiMiC Communication Library, enables communication between external programs
    coupled through the MiMiC framework. See https://mimic-project.org/ for further information.
    """

    homepage = "https://gitlab.com/mimic-project/mcl"
    url = "https://gitlab.com/mimic-project/mcl/-/archive/3.0.0/mcl-3.0.0.tar.gz"
    git = "https://gitlab.com/mimic-project/mcl.git"

    license("GPL-3.0", checked_by="mkrack")

    maintainers("mkrack")

    version("3.0.0", sha256="3e740582836fe90e04a693cfc5a219826bcac03217f70ea5570bad6aeafda685")

    variant("build_tests", default=False, description="Build tests")
    variant("coverage", default=False, description="Enable code coverage report")
    variant("build_fortran_api", default=True, description="Build Fortran API module")
    variant("build_shared_libs", default=True, description="Build using shared libraries")
    variant(
        "disable_mpi_f08",
        default=False,
        description="Disable MPI F08 Fortran module (even if it is available)",
    )
    variant(
        "build_type",
        default="Release",
        description="CMake build type",
        values=("Debug", "Release"),
    )

    depends_on("c", type="build")
    depends_on("cxx", type="build")

    depends_on("fortran", when="+build_fortran_api", type="build")

    depends_on("cmake@3.12:", type="build")
    depends_on("mpi@2.1:", type="build")
    depends_on("python@3.6:", type="build")

    build_targets = ["mcl"]

    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_TESTS", "build_tests"),
            self.define_from_variant("ENABLE_COVERAGE", "coverage"),
            self.define_from_variant("BUILD_FORTRAN_API", "build_fortran_api"),
            self.define_from_variant("BUILD_SHARED_LIBS", "build_shared_libs"),
            self.define_from_variant("DISABLE_MPI_F08", "disable_mpi_f08"),
        ]
        return args
