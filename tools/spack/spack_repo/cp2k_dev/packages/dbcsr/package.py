# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os

from spack_repo.builtin.build_systems.cmake import CMakePackage, generator
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage

from spack.package import *


class Dbcsr(CMakePackage, CudaPackage, ROCmPackage):
    """Distributed Block Compressed Sparse Row matrix library."""

    homepage = "https://github.com/cp2k/dbcsr"
    git = "https://github.com/cp2k/dbcsr.git"
    url = "https://github.com/cp2k/dbcsr/releases/download/v2.2.0/dbcsr-2.2.0.tar.gz"
    list_url = "https://github.com/cp2k/dbcsr/releases"

    maintainers("dev-zero", "mtaillefumier", "RMeli", "hfp")

    license("GPL-2.0-or-later")

    version("develop", branch="develop")
    version("2.9.2-preview", commit="4d85b72e3427ec7f595ad0cf38a4d00f72f9ff60")
    version("2.9.1", sha256="fa5a4aeba0a07761511af2c26c779bd811b5ea0ef06a5d94535b6dd7b2e0ce59")
    version("2.9.0", sha256="a04cacd2203bd97a31ac993f9ab84237a48191140bba29efadbc27db544bbcd6")
    version("2.8.0", sha256="d55e4f052f28d1ed0faeaa07557241439243287a184d1fd27f875c8b9ca6bd96")

    new_smm_versions = "@2.9.2-preview,2.9.2:,develop"

    variant("tests", default=False, description="Build DBCSR unit tests")
    variant("tests", default=True, description="Build DBCSR unit tests", when="@2.1:2.2")
    variant("mpi", default=True, description="Compile with MPI")
    variant("openmp", default=True, description="Build with OpenMP support")
    variant("shared", default=True, description="Build shared library")
    variant(
        "smm",
        default="libxsmm",
        values=("libxsmm", "libxs", "blas"),
        description="Library for small matrix multiplications",
    )
    variant(
        "cuda_arch_35_k20x",
        default=False,
        description=(
            "CP2K (resp. DBCSR) has specific parameter sets for"
            " different GPU models. Enable this when building"
            " with cuda_arch=35 for a K20x instead of a K40"
        ),
    )
    variant("examples", default=True, description="Build examples")

    variant("opencl", default=False, description="Enable OpenCL backend")
    variant("mpi_f08", default=False, when="@2.6:", description="Use mpi F08 module")

    depends_on("c", type="build")  # generated
    depends_on("cxx", type="build")  # generated
    depends_on("fortran", type="build")  # generated

    depends_on("blas")
    depends_on("lapack")
    depends_on("mpi", when="+mpi")

    conflicts(
        "smm=libxs", when="@:2.9.1", msg="smm=libxs is only supported by DBCSR 2.9.2 or newer"
    )

    depends_on("libxsmm@1.11:", when="@:2.9.1 smm=libxsmm")

    with when(f"{new_smm_versions} smm=libxs"):
        depends_on("libxs")

    with when(f"{new_smm_versions} smm=libxsmm"):
        depends_on("libxs")
        depends_on("libxsmm")

    depends_on("cmake@3.10:", type="build")
    depends_on("cmake@3.12:", type="build", when="@2.1:")
    depends_on("cmake@3.17:", type="build", when="@2.2:")
    depends_on("cmake@3.22:", type="build", when="@2.3:")

    depends_on("py-fypp", type="build")
    depends_on("py-fypp@3.1:", type="build", when="@2.6:")
    depends_on("pkgconfig", type="build")
    depends_on("python@3.6:", type="build", when="+cuda")

    depends_on("hipblas", when="+rocm")

    # Several packages provide "opencl" (incl. ICD/loader), e.g., "cuda"
    depends_on("opencl", when="+opencl")
    depends_on("libxstream", when="+opencl")
    opencl_loader_header_version = "2022.10.24"
    depends_on(f"opencl-c-headers@{opencl_loader_header_version}:", when="+opencl")
    requires(f"%opencl=opencl-icd-loader@{opencl_loader_header_version}:", when="+opencl")

    # All examples require MPI
    conflicts("+examples", when="~mpi", msg="Examples require MPI")

    # We only support specific gpu archs for which we have parameter files
    # for optimal kernels. Note that we don't override the parent class arch
    # properties, since the parent class defines constraints for different archs
    # Instead just mark all unsupported cuda archs as conflicting.
    dbcsr_cuda_archs = ("35", "37", "60", "70", "80", "90")
    cuda_msg = f"dbcsr only supports cuda_arch {dbcsr_cuda_archs}"

    for arch in CudaPackage.cuda_arch_values:
        if arch not in dbcsr_cuda_archs:
            conflicts("+cuda", when=f"cuda_arch={arch}", msg=cuda_msg)

    conflicts("+cuda", when="cuda_arch=none", msg=cuda_msg)

    dbcsr_amdgpu_targets = (
        "gfx906",
        "gfx908",
        "gfx90a",
        "gfx90a:xnack-",
        "gfx90a:xnack+",
        "gfx942",
        "gfx950",
    )
    amd_msg = f"""DBCSR supports these AMD gpu targets:  {", ".join(dbcsr_amdgpu_targets)}.
                  Set amdgpu_target explicitly to one of the supported targets"""

    for arch in ROCmPackage.amdgpu_targets:
        if arch not in dbcsr_amdgpu_targets:
            conflicts("+rocm", when=f"amdgpu_target={arch}", msg=amd_msg)

    # GPU runtimes are usually mutually exclusive
    accel_msg = "CUDA and ROCm are mutually exlusive"
    conflicts("+cuda", when="+rocm", msg=accel_msg)
    conflicts("+cuda", when="+opencl", msg=accel_msg)
    conflicts("+rocm", when="+opencl", msg=accel_msg)

    # Require OpenMP threading by making other options conflict
    conflicts("^intel-oneapi-mkl threads=none", when="+openmp")
    conflicts("^intel-oneapi-mkl threads=tbb", when="+openmp")
    conflicts("^openblas threads=pthreads", when="+openmp")
    conflicts("^openblas threads=none", when="+openmp")

    # OpenCL follows the legacy LIBXSMM SMM path up to DBCSR 2.9.1.
    requires("smm=libxsmm", when="+opencl @:2.9.1")

    # Starting with the 2.9.2 interface, LIBXS is an explicit SMM backend.
    requires("smm=libxs", when=f"+opencl {new_smm_versions}")

    with when("+mpi"):
        # When using mpich 4.1 or higher, mpi_f08 has to be used, otherwise:
        # Error: Type mismatch in argument 'baseptr' at (1); passed TYPE(c_ptr)
        # to INTEGER(8)
        conflicts("^mpich@4.1:", when="@:2.5")
        conflicts("~mpi_f08", when="^mpich@4.1:")
        depends_on("mpich+fortran", when="^[virtuals=mpi] mpich")

    generator("ninja")
    depends_on("ninja@1.10:", type="build")

    @when("+rocm")
    def patch(self):
        for directory, subdirectory, files in os.walk(os.getcwd()):
            for i in files:
                file_path = os.path.join(directory, i)
                filter_file("USE ISO_C_BINDING", "USE,INTRINSIC :: ISO_C_BINDING", file_path)
                filter_file("USE ISO_FORTRAN_ENV", "USE,INTRINSIC :: ISO_FORTRAN_ENV", file_path)
                filter_file("USE omp_lib", "USE,INTRINSIC :: omp_lib", file_path)
                filter_file("USE OMP_LIB", "USE,INTRINSIC :: OMP_LIB", file_path)
                filter_file("USE iso_c_binding", "USE,INTRINSIC :: iso_c_binding", file_path)
                filter_file("USE iso_fortran_env", "USE,INTRINSIC :: iso_fortran_env", file_path)

    def cmake_args(self):
        spec = self.spec

        if "+cuda" in spec and len(spec.variants["cuda_arch"].value) > 1:
            raise InstallError("dbcsr supports only one cuda_arch at a time")

        if "+rocm" in spec and len(spec.variants["amdgpu_target"].value) > 1:
            raise InstallError("DBCSR supports only one amdgpu_arch at a time")

        smm = spec.variants["smm"].value
        has_libxs_smm = (
            spec.satisfies("@2.9.2-preview")
            or spec.satisfies("@2.9.2:")
            or spec.satisfies("@develop")
        )

        args = [
            self.define_from_variant("USE_MPI", "mpi"),
            self.define_from_variant("USE_OPENMP", "openmp"),
            # C API needs MPI
            self.define_from_variant("WITH_C_API", "mpi"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("WITH_EXAMPLES", "examples"),
            self.define_from_variant("BUILD_TESTING", "tests"),
        ]

        if has_libxs_smm:
            args += [
                self.define("USE_LIBXS", smm in ("libxs", "libxsmm")),
                self.define("USE_LIBXSMM", smm == "libxsmm"),
            ]
        else:
            args.append("-DUSE_SMM=%s" % ("libxsmm" if smm == "libxsmm" else "blas"))

        lapack, blas = spec["lapack"], spec["blas"]
        if blas.name != "intel-oneapi-mkl":
            args += [
                "-DBLAS_FOUND=true",
                "-DBLAS_LIBRARIES=%s" % (blas.libs.joined(";")),
                "-DLAPACK_FOUND=true",
                "-DLAPACK_LIBRARIES=%s" % (lapack.libs.joined(";")),
            ]

        if self.spec.satisfies("+cuda"):
            cuda_arch = self.spec.variants["cuda_arch"].value[0]

            gpu_map = {
                "35": "K40",
                "37": "K80",
                "60": "P100",
                "70": "V100",
                "80": "A100",
                "90": "H100",
            }

            gpuver = gpu_map[cuda_arch]
            if cuda_arch == "35" and self.spec.satisfies("+cuda_arch_35_k20x"):
                gpuver = "K20X"

            args += ["-DWITH_GPU=%s" % gpuver, "-DUSE_ACCEL=cuda"]

        if self.spec.satisfies("+rocm"):
            amd_arch = self.spec.variants["amdgpu_target"].value[0]
            gpuver = {
                "gfx906": "Mi50",
                "gfx908": "Mi100",
                "gfx90a": "Mi250",
                "gfx90a:xnack-": "Mi250",
                "gfx90a:xnack+": "Mi250",
                "gfx942": "Mi300",
                "gfx950": "Mi350",
            }[amd_arch]

            args += [f"-DWITH_GPU={gpuver}", "-DUSE_ACCEL=hip"]

        if self.spec.satisfies("+opencl"):
            args += ["-DUSE_ACCEL=opencl"]

        if self.spec.satisfies("+mpi_f08"):
            args += ["-DUSE_MPI_F08=ON"]

        return args

    def check(self):
        """Override CMakePackage's check() to enforce seralized test runs
        since they are already parallelized"""
        with working_dir(self.build_directory):
            self._if_ninja_target_execute("test", parallel=False)
