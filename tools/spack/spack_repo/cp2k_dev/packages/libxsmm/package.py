# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libxsmm(CMakePackage):
    """Library for specialized dense
    and sparse matrix operations,
    and deep learning primitives."""

    homepage = "https://github.com/libxsmm/libxsmm"
    url = "https://github.com/libxsmm/libxsmm/archive/1.17.tar.gz"
    git = "https://github.com/libxsmm/libxsmm.git"

    maintainers("hfp")

    license("BSD-3-Clause")

    version("20260617", commit="db07b74b2ffcf9f303bc0ad0b14740512c5341d1")
    version("20260610-cp2k", commit="79033a7a1da039681cfa3b5ed2fe662a401fa4f3")
    version("20260526-cp2k", commit="0cea22fdc34ec54bc59ffb47a43cb3e28b26d3e0")
    version("1.17-cp2k", commit="e0c4a2389afba36c453233ad7de07bd92c715bec")

    variant("shared", default=False, description="With shared libraries.")
    variant("fortran", default=True, description="With Fortran support.")

    depends_on("c", type="build")  # generated
    depends_on("cxx", type="build")  # generated
    depends_on("fortran", type="build", when="+fortran")  # generated

    depends_on("python", type="build")
    depends_on("cmake@3.13:", type="build")

    # A recent `as` is needed to compile libxmss until version 1.17
    # (<https://github.com/spack/spack/issues/28404>), but not afterwards
    # (<https://github.com/spack/spack/pull/21671#issuecomment-779882282>).
    # depends_on("binutils+ld+gas@2.33:", type="build")

    # Version 2.0 supports both x86_64 and aarch64
    requires("target=x86_64:", "target=aarch64:")
    # requires("target=x86_64:", when="@:1")

    @property
    def libs(self):
        result = find_libraries(["libxsmm", "libxsmmf"], root=self.prefix, recursive=True)
        if len(result) == 0:
            result = find_libraries(
                ["libxsmm", "libxsmmf"], root=self.prefix, shared=False, recursive=True
            )
        return result

    def cmake_args(self):
        spec = self.spec
        return [
            self.define("BUILD_SHARED_LIBS", spec.satisfies("+shared")),
            self.define("LIBXSMM_FORTRAN", spec.satisfies("+fortran")),
        ]
