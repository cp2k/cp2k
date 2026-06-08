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

    maintainers("awvwgk", "mkrack")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("dev20260522", branch="skala", commit="69ee054054c642bce0d9b5e8c9c0c13afa97b774")
    version(
        "1.2.dev2",
        sha256="6659d00522b443f2557a960a9ba0217449d2ee25df809634d6becca23eb0f1ff",
        url="https://github.com/microsoft/skala/releases/download/v1.1.1/gauxc-skala-r2.tar.gz",
    )
    version("1.1", sha256="17de077fb23e44d03b0ed14dcd8117c01e5b3431fbefa2352d751639ada7f91c")

    variant("c", default=True, description="Build with C API support")
    variant("fortran", default=True, description="Build with Fortran support")
    variant("host", default=True, description="Build with host integrator")
    variant("cuda", default=False, description="Build with CUDA support")
    variant("hip", default=False, description="Build with HIP support")
    variant("mpi", default=False, description="Build with MPI support")
    variant("openmp", default=True, description="Build with OpenMP support")
    variant("tests", default=False, description="Build with unit tests")
    variant("gau2grid", default=True, description="Build with Gau2Grid collocation")
    variant("hdf5", default=True, description="Build with HDF5 support")
    variant("onedft", default=False, description="Build with ONEDFT functional")
    variant("skala", default=False, description="Build with Skala support")
    variant("fast_rsqrt", default=False, description="Build with fast RSQRT support")
    variant("blas_prefer_ilp64", default=False, description="Prefer ILP64 for host BLAS")
    variant("link_cuda_static", default=False, description="Link GauXC with static CUDA libs")
    variant("magma", default=False, description="Build with MAGMA linear algebra support")
    variant("nccl", default=False, description="Build with NCCL collectives")
    variant("cutlass", default=False, description="Build with CUTLASS linear algebra support")
    variant("pic", default=True, description="Build position independent code")

    conflicts("+skala", when="~onedft", msg="Skala requires ONEDFT support")
    conflicts("+magma", when="~cuda ~hip", msg="MAGMA requires CUDA or HIP")
    conflicts("+nccl", when="~cuda", msg="NCCL requires CUDA")
    conflicts("+cutlass", when="~cuda", msg="CUTLASS requires CUDA")
    conflicts("+fortran", when="~c", msg="Fortran bindings require C API")

    # Skala resource file
    skala_version = "1.1"
    skala_file = f"skala-{skala_version}.fun"
    resource(
        name="skala_fun",
        url=f"https://huggingface.co/microsoft/skala-{skala_version}/resolve/main/{skala_file}",
        sha256="0c8432ac3f03c8f1276372df9aca5b7ee7f8939d47a8789eb158976e89aa0606",
        expand=False,
        placement=skala_file,
        when="+skala",
    )

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", when="+fortran", type="build")
    depends_on("cmake@3.20:", type="build")
    depends_on("ninja@1.10:", type="build")
    depends_on("exchcxx~cuda", when="~cuda")
    depends_on("exchcxx+cuda", when="+cuda")
    depends_on("integratorxx")
    depends_on("blas")
    depends_on("hdf5", when="+hdf5")
    depends_on("highfive@:2.10.1", when="+hdf5")
    depends_on("mpi", when="+mpi")
    depends_on("cuda", when="+cuda")
    depends_on("hip", when="+hip")
    depends_on("gau2grid", when="+gau2grid")
    depends_on("nlohmann-json", when="+skala")
    depends_on("py-torch@2:", when="+skala")
    depends_on("py-torch@2:+cuda", when="+skala+cuda")
    depends_on("magma", when="+magma")
    depends_on("nccl", when="+nccl")
    depends_on("cutlass", when="+cutlass")

    generator("ninja")

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define("GAUXC_ENABLE_C", spec.satisfies("+c")),
            self.define("GAUXC_ENABLE_FORTRAN", spec.satisfies("+fortran")),
            self.define("GAUXC_ENABLE_HOST", spec.satisfies("+host")),
            self.define("GAUXC_ENABLE_CUDA", spec.satisfies("+cuda")),
            self.define("GAUXC_ENABLE_HIP", spec.satisfies("+hip")),
            self.define("GAUXC_ENABLE_MPI", spec.satisfies("+mpi")),
            self.define("GAUXC_ENABLE_OPENMP", spec.satisfies("+openmp")),
            self.define("GAUXC_ENABLE_TESTS", spec.satisfies("+tests")),
            self.define("GAUXC_ENABLE_GAU2GRID", spec.satisfies("+gau2grid")),
            self.define("GAUXC_ENABLE_HDF5", spec.satisfies("+hdf5")),
            self.define("GAUXC_ENABLE_ONEDFT", spec.satisfies("+onedft")),
            self.define("GAUXC_ENABLE_FAST_RSQRT", spec.satisfies("+fast_rsqrt")),
            self.define("GAUXC_BLAS_PREFER_ILP64", spec.satisfies("+blas_prefer_ilp64")),
            self.define("GAUXC_LINK_CUDA_STATIC", spec.satisfies("+link_cuda_static")),
            self.define("GAUXC_ENABLE_MAGMA", spec.satisfies("+magma")),
            self.define("GAUXC_ENABLE_NCCL", spec.satisfies("+nccl")),
            self.define("GAUXC_ENABLE_CUTLASS", spec.satisfies("+cutlass")),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        return args

    def install(self, spec, prefix):
        super().install(spec, prefix)
        if self.spec.satisfies("+skala"):
            target_dir = self.prefix.share.gauxc.join("onedft_models")
            mkdirp(target_dir)
            print(target_dir)
            install(join_path(self.skala_file, self.skala_file), target_dir)
