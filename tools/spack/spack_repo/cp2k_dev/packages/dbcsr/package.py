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

    # homepage = "https://github.com/cp2k/dbcsr"
    # git = "https://github.com/cp2k/dbcsr.git"
    # url = "https://github.com/cp2k/dbcsr/releases/download/v2.2.0/dbcsr-2.2.0.tar.gz"
    # list_url = "https://github.com/cp2k/dbcsr/releases"
    homepage = "https://github.com/Growl1234/dbcsr"
    git = "https://github.com/Growl1234/dbcsr.git"

    maintainers("dev-zero", "mtaillefumier", "RMeli", "hfp")

    license("GPL-2.0-or-later")

    version("develop", branch="develop")
    version("cmake-test", commit="4ec13670d719f580e7fe73cd1d263ad9558ced66")

    variant("tests", default=False, description="Build DBCSR unit tests")
    variant("tests", default=True, description="Build DBCSR unit tests", when="@2.1:2.2")
    variant("mpi", default=True, description="Compile with MPI")
    variant("openmp", default=True, description="Build with OpenMP support")
    variant("shared", default=True, description="Build shared library")
    variant(
        "smm",
        default="libxsmm",
        values=("libxsmm", "blas"),
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

    variant("g2g", default=False, when="@:2.8", description="GPU-aware MPI with CUDA/HIP")
    conflicts("+g2g", when="~cuda ~rocm", msg="GPU-aware MPI requires +cuda or +rocm")

    depends_on("c", type="build")  # generated
    depends_on("cxx", type="build")  # generated
    depends_on("fortran", type="build")  # generated

    depends_on("blas")
    depends_on("lapack")
    depends_on("mpi", when="+mpi")

    with when("smm=libxsmm"):
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
        "gfx910",
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

    # OpenCL backend implementation relies on LIBXSMM
    requires("smm=libxsmm", when="+opencl")

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

        args = [
            "-DUSE_SMM=%s" % ("libxsmm" if "smm=libxsmm" in spec else "blas"),
            self.define_from_variant("USE_MPI", "mpi"),
            self.define_from_variant("USE_OPENMP", "openmp"),
            # C API needs MPI
            self.define_from_variant("WITH_C_API", "mpi"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("WITH_EXAMPLES", "examples"),
            self.define_from_variant("WITH_G2G", "g2g"),
            self.define_from_variant("BUILD_TESTING", "tests"),
        ]

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
