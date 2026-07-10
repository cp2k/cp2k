# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage

from spack.package import *


class Gauxc(CMakePackage, CudaPackage, ROCmPackage):
    """GauXC is a modern, modular C++ library for the evaluation of quantities related to
    the exchange-correlation (XC) energy (e.g. potential, etc) in the Gaussian basis set
    discretization of Kohn-Sham density function theory (KS-DFT) on heterogenous
    architectures."""

    homepage = "https://github.com/wavefunction91/GauXC"
    git = "https://github.com/wavefunction91/GauXC.git"
    url = "https://github.com/wavefunction91/GauXC/archive/refs/tags/v0.0.tar.gz"

    maintainers("mkrack", "RMeli")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("dev20260608", branch="skala", commit="36aeab80bf7b2245087bdfeffd8e61b11206b3bd")
    version(
        "1.2.dev2",
        sha256="6659d00522b443f2557a960a9ba0217449d2ee25df809634d6becca23eb0f1ff",
        url="https://github.com/microsoft/skala/releases/download/v1.1.1/gauxc-skala-r2.tar.gz",
    )
    version("1.1", sha256="17de077fb23e44d03b0ed14dcd8117c01e5b3431fbefa2352d751639ada7f91c")

    variant("c", default=True, description="Build with C API support", when="@master")
    variant("c", default=True, description="Build with C API support", when="@1.2.dev2")
    variant("c", default=True, description="Build with C API support", when="@dev20260608")
    variant("fortran", default=True, description="Build with Fortran support", when="@dev20260608")
    variant("host", default=True, description="Build with host integrator")
    variant("mpi", default=False, description="Build with MPI support")
    variant("openmp", default=True, description="Build with OpenMP support")
    variant("tests", default=False, description="Build with unit tests")
    variant("hdf5", default=True, description="Build with HDF5 support")
    variant("skala", default=False, description="Build with Skala support")
    variant("fast_rsqrt", default=False, description="Build with fast RSQRT support")
    variant("blas_prefer_ilp64", default=False, description="Prefer ILP64 for host BLAS")
    variant("link_cuda_static", default=False, description="Link GauXC with static CUDA libs")
    variant("magma", default=False, description="Build with MAGMA linear algebra support")
    variant("nccl", default=False, description="Build with NCCL collectives")
    variant("cutlass", default=False, description="Build with CUTLASS linear algebra support")
    variant("pic", default=True, description="Build position independent code")

    # Skala resource file
    SKALA_RESOURCES = {"1.1": "0c8432ac3f03c8f1276372df9aca5b7ee7f8939d47a8789eb158976e89aa0606"}
    variant("skala_version", default="none", values=("none", *SKALA_RESOURCES.keys()), multi=False)
    for skala_version, sha256 in SKALA_RESOURCES.items():
        resource(
            name=f"skala_fun_{skala_version}",
            url=f"https://huggingface.co/microsoft/skala-{skala_version}/resolve/main/skala-{skala_version}.fun",
            sha256=sha256,
            expand=False,
            placement="onedft_models",
            when=f"+skala skala_version={skala_version}",
        )

    conflicts("+magma", when="~cuda ~rocm", msg="MAGMA requires CUDA or HIP")
    conflicts("+nccl", when="~cuda", msg="NCCL requires CUDA")
    conflicts("+cutlass", when="~cuda", msg="CUTLASS requires CUDA")
    conflicts("+fortran", when="~c @1.2.dev2:", msg="Fortran bindings require C API")

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
    depends_on("hip", when="+rocm")
    depends_on("gau2grid")
    depends_on("nlohmann-json", when="+skala")
    depends_on("py-torch@2:~cuda", when="+skala~cuda")
    depends_on("py-torch@2:+cuda", when="+skala+cuda")
    depends_on("magma", when="+magma")
    depends_on("nccl", when="+nccl")
    depends_on("cutlass", when="+cutlass")

    generator("ninja")

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define_from_variant("GAUXC_ENABLE_C", "c"),
            self.define_from_variant("GAUXC_ENABLE_FORTRAN", "fortran"),
            self.define_from_variant("GAUXC_ENABLE_HOST", "host"),
            self.define_from_variant("GAUXC_ENABLE_CUDA", "cuda"),
            self.define_from_variant("GAUXC_ENABLE_HIP", "rocm"),
            self.define_from_variant("GAUXC_ENABLE_MPI", "mpi"),
            self.define_from_variant("GAUXC_ENABLE_OPENMP", "openmp"),
            self.define_from_variant("GAUXC_ENABLE_TESTS", "tests"),
            self.define_from_variant("GAUXC_ENABLE_HDF5", "hdf5"),
            self.define_from_variant("GAUXC_ENABLE_ONEDFT", "skala"),
            self.define_from_variant("GAUXC_ENABLE_FAST_RSQRT", "fast_rsqrt"),
            self.define_from_variant("GAUXC_BLAS_PREFER_ILP64", "blas_prefer_ilp64"),
            self.define_from_variant("GAUXC_LINK_CUDA_STATIC", "link_cuda_static"),
            self.define_from_variant("GAUXC_ENABLE_MAGMA", "magma"),
            self.define_from_variant("GAUXC_ENABLE_NCCL", "nccl"),
            self.define_from_variant("GAUXC_ENABLE_CUTLASS", "cutlass"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        if spec.satisfies("+cuda"):
            archs = spec.variants["cuda_arch"].value
            if "none" not in archs:
                arch_str = ";".join(archs)
                args.append(self.define("CMAKE_CUDA_ARCHITECTURES", arch_str))
        if spec.satisfies("+rocm"):
            archs = spec.variants["amdgpu_target"].value
            if "none" not in archs:
                arch_str = ";".join(archs)
                args.append(self.define("CMAKE_HIP_ARCHITECTURES", arch_str))
        return args

    def install(self, spec, prefix):
        super().install(spec, prefix)
        if spec.satisfies("+skala"):
            skala_version = spec.variants["skala_version"].value
            skala_file = f"skala-{skala_version}.fun"
            target_dir = prefix.share.gauxc.join("onedft_models")
            mkdirp(target_dir)
            install(join_path("onedft_models", skala_file), target_dir)
