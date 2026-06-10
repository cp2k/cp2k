# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libxsmm(CMakePackage):
    """Library for specialized dense
    and sparse matrix operations,
    and deep learning primitives."""

    # homepage = "https://github.com/libxsmm/libxsmm"
    # url = "https://github.com/libxsmm/libxsmm/archive/1.17.tar.gz"
    # git = "https://github.com/libxsmm/libxsmm.git"
    homepage = "https://github.com/Growl1234/libxsmm"
    url = "https://github.com/Growl1234/libxsmm/archive/1.17.tar.gz"
    git = "https://github.com/Growl1234/libxsmm.git"

    maintainers("hfp")

    license("BSD-3-Clause")

    version("cmake-test", commit="b35d767abffc7bca845a81f7e408fec68300821c")
    version("20260526-cp2k", commit="0cea22fdc34ec54bc59ffb47a43cb3e28b26d3e0")
    version("1.17-cp2k", commit="e0c4a2389afba36c453233ad7de07bd92c715bec")

    variant("shared", default=False, description="With shared libraries (and static libraries).")
    variant("debug", default=False, description="With call-trace (LIBXSMM_TRACE); unoptimized.")
    variant(
        "blas",
        default="default",
        multi=False,
        description="Control behavior of BLAS calls",
        values=("default", "0", "1", "2"),
    )

    depends_on("c", type="build")  # generated
    depends_on("cxx", type="build")  # generated
    depends_on("fortran", type="build")  # generated

    depends_on("python", type="build")
    depends_on("cmake@3.13:", type="build")
    depends_on("blas", when="blas=1")
    depends_on("blas", when="blas=2")

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
        blas_val = spec.variants["blas"].value

        return [
            self.define("XSMM_STATIC", not spec.satisfies("+shared")),
            self.define("XSMM_BLAS", blas_val in ("1", "2")),
        ]
