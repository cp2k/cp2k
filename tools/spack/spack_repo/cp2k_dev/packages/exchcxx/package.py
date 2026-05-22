# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack_repo.builtin.build_systems.cuda import CudaPackage

from spack.package import *


class Exchcxx(CMakePackage, CudaPackage):
    """Exchange correlation (XC) library for density functional theory (DFT)
    calculations in modern C++."""

    homepage = "https://github.com/wavefunction91/ExchCXX"
    git = "https://github.com/wavefunction91/ExchCXX.git"
    url = "https://github.com/wavefunction91/ExchCXX/archive/refs/tags/v1.0.0.tar.gz"

    maintainers("awvwgk")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("1.0.0", sha256="0ea38aab563ea5d11a03d80dc1bf909821091d903eb852c5cb350484753a847a")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("cmake@3.20:", type="build")
    depends_on("ninja@1.10:", type="build")
    depends_on("libxc@7.0.0:")
    depends_on("cuda", when="+cuda")

    generator("ninja")

    variant("cuda", default=False, description="Build with CUDA support")

    def cmake_args(self):
        args = []
        if "+cuda" in self.spec:
            args.append(self.define("EXCHCXX_ENABLE_CUDA", True))
        return args
