# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack_repo.builtin.build_systems.cuda import CudaPackage

from spack.package import *


class Gauxc(CMakePackage, CudaPackage):
    """GauXC is a modern, modular C++ library for the evaluation of quantities related to
    the exchange-correlation (XC) energy (e.g. potential, etc) in the Gaussian basis set
    discretization of Kohn-Sham density function theory (KS-DFT) on heterogenous architectures."""

    homepage = "https://github.com/wavefunction91/GauXC"
    git = "https://github.com/wavefunction91/GauXC.git"
    url = "https://github.com/wavefunction91/GauXC/archive/refs/tags/v1.1.tar.gz"

    maintainers("awvwgk")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("dev20260522", branch="skala", commit="69ee054054c642bce0d9b5e8c9c0c13afa97b774")
    version(
        "1.2.dev2",
        sha256="6659d00522b443f2557a960a9ba0217449d2ee25df809634d6becca23eb0f1ff",
        url="https://github.com/microsoft/skala/releases/download/v1.1.1/gauxc-skala-r2.tar.gz",
    )
    version("1.1", sha256="17de077fb23e44d03b0ed14dcd8117c01e5b3431fbefa2352d751639ada7f91c")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", when="+fortran", type="build")
    depends_on("cmake@3.20:", type="build")
    depends_on("ninja@1.10:", type="build")
    depends_on("exchcxx ~cuda", when="~cuda")
    depends_on("exchcxx +cuda", when="+cuda")
    depends_on("integratorxx")
    depends_on("hdf5", when="+hdf5")
    depends_on("mpi", when="+mpi")
    depends_on("cuda", when="+cuda")
    depends_on("highfive@:2.10.1", when="+hdf5")
    depends_on("gau2grid")
    depends_on("blas")
    depends_on("nlohmann-json", when="+skala")
    depends_on("py-torch@2:", when="+skala")
    depends_on("py-torch@2: +cuda", when="+skala+cuda")

    generator("ninja")

    variant("openmp", default=True, description="Build with OpenMP support")
    variant("mpi", default=False, description="Build with MPI support")
    variant("cuda", default=False, description="Build with CUDA support")
    variant("hdf5", default=False, description="Build with HDF5 support")
    variant("c", default=False, description="Build with C support")
    variant("fortran", default=False, description="Build with Fortran support")
    variant("skala", default=False, description="Build with Skala support")
    variant("host", default=True, description="Build with host integrators")
    variant("tests", default=False, description="Build with testing framework (Catch2)")
    variant("pic", default=True, description="Build position independent code")

    def cmake_args(self):
        args = [
            self.define("GAUXC_ENABLE_OPENMP", self.spec.satisfies("+openmp")),
            self.define("GAUXC_ENABLE_MPI", self.spec.satisfies("+mpi")),
            self.define("GAUXC_ENABLE_CUDA", self.spec.satisfies("+cuda")),
            self.define("GAUXC_ENABLE_HDF5", self.spec.satisfies("+hdf5")),
            self.define("GAUXC_ENABLE_C", self.spec.satisfies("+c")),
            self.define("GAUXC_ENABLE_FORTRAN", self.spec.satisfies("+fortran")),
            self.define("GAUXC_ENABLE_ONEDFT", self.spec.satisfies("+skala")),
            self.define("GAUXC_ENABLE_HOST", self.spec.satisfies("+host")),
            self.define("GAUXC_ENABLE_TESTS", self.spec.satisfies("+tests")),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        return args
