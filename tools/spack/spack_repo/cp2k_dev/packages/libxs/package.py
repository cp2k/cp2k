# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libxs(CMakePackage):
    """LIBXS is a portable C library providing building blocks for memory
    operations, numerics, synchronization, and more -- with a focus on
    performance and minimal dependencies. Targets x86-64, AArch64, and RISC-V;
    requires only a C89 compiler. Originally developed as part of LIBXSMM.
    """

    maintainers("hfp", "mkrack")

    homepage = "hhttps://libxs.readthedocs.io/en/latest"
    url = (
        "https://github.com/hfp/libxs/archiv/libxs-ab416130f8c9f7edb8c1bf3d3abaf402f61d0fe0.tar.gz"
    )
    git = "https://github.com/hfp/libxs.git"

    license("BSD-3-Clause", checked_by="mkrack")

    version("main", branch="main")
    version("20260605", commit="ab416130f8c9f7edb8c1bf3d3abaf402f61d0fe0")

    variant("fortran", default=True, description="Build Fortran module interface")
    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=True, description="Build shared libraries (otherwise static)")

    depends_on("cmake@3.13:", type="build")
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build")

    depends_on("libxsmm")

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define_from_variant("LIBXS_FORTRAN", "fortran"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define_from_variant("LIBXS_SHARED", "shared"),
        ]
        return args
