# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libxs(CMakePackage):
    """LIBXS is a portable C library providing building blocks for memory
    operations, numerics, synchronization, and more -- with a focus on
    performance and minimal dependencies. Targets x86-64, AArch64, and RISC-V;
    requires only a C89 compiler. Originally developed as part of LIBXSMM."""

    homepage = "https://libxs.readthedocs.io/en/latest"
    git = "https://github.com/hfp/libxs.git"
    url = "https://github.com/hfp/libxs/releases/download/1.0.0/libxs-1.0.0.tar.gz"

    maintainers("hfp", "mkrack", "mtaillefumier")

    license("BSD-3-Clause", checked_by="mkrack")

    version("main", branch="main")
    version("1.0.0", sha256="de26f50cb986a2f0e4f92c0eb489d40a44f7e4c5acd22751a6cfa2829dabd04d")

    variant("fortran", default=False, description="Build Fortran module interface")
    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=False, description="Build shared libraries (otherwise static)")
    variant("tests", default=False, description="Build unit testsi (requires BLAS)")

    depends_on("cmake@3.13:", type="build")
    depends_on("c", type="build")
    depends_on("fortran", type="build", when="+fortran")

    depends_on("blas", when="+tests")

    def cmake_args(self):
        args = [
            self.define_from_variant("LIBXS_FORTRAN", "fortran"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("BUILD_TESTING", "tests"),
        ]
        return args
