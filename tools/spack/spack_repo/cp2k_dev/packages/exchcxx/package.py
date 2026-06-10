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

    maintainers("awvwgk", "mkrack")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("1.0.0", sha256="0ea38aab563ea5d11a03d80dc1bf909821091d903eb852c5cb350484753a847a")

    variant("benchmark", default=False, description="Enable benchmarking")
    variant("cuda", default=False, description="Build with CUDA support")
    variant("hip", default=False, description="Build with HIP support")
    variant("sycl", default=False, description="Build with SYCL support")
    variant("libxc", default=True, description="Build with LibXC support")
    variant("shared", default=False, description="Build shared libraries")
    variant("pic", default=True, description="Build position independent code")

    conflicts("+cuda", when="+hip", msg="CUDA and HIP are mutually exclusive")
    conflicts("+cuda", when="+sycl", msg="CUDA and SYCL are mutually exclusive")
    conflicts("+hip", when="+sycl", msg="HIP and SYCL are mutually exclusive")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("cmake@3.21:", type="build")
    depends_on("ninja@1.10:", type="build")
    depends_on("libxc@7.0.0:", when="+libxc")
    depends_on("cuda@11:", when="+cuda")
    depends_on("hip", when="+hip")
    depends_on("intel-oneapi-dpcpp", when="+sycl")

    generator("ninja")

    def cmake_args(self):
        args = [
            self.define_from_variant("EXCHCXX_ENABLE_BENCHMARK", "benchmark"),
            self.define_from_variant("EXCHCXX_ENABLE_CUDA", "cuda"),
            self.define_from_variant("EXCHCXX_ENABLE_HIP", "hip"),
            self.define_from_variant("EXCHCXX_ENABLE_SYCL", "sycl"),
            self.define_from_variant("EXCHCXX_ENABLE_LIBXC", "libxc"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        if self.spec.satisfies("+cuda"):
            args.append(CMakeBuilder.define_cuda_architectures(self))
        return args
