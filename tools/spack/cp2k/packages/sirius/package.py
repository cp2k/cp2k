# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os

from spack.package import *


class Sirius(CMakePackage, CudaPackage, ROCmPackage):
    """Domain specific library for electronic structure calculations"""

    homepage = "https://github.com/electronic-structure/SIRIUS"
    url = "https://github.com/electronic-structure/SIRIUS/archive/v6.1.5.tar.gz"
    list_url = "https://github.com/electronic-structure/SIRIUS/releases"
    git = "https://github.com/electronic-structure/SIRIUS.git"

    maintainers("simonpintarelli", "haampie", "dev-zero", "AdhocMan", "toxa81", "RMeli")

    license("BSD-2-Clause")

    version("develop", branch="develop")
    version("master", branch="master")

    version(
        "7.7.0",
        sha256="be0bdc76db9eb8afdcb950f0ccaf7535b8e85d72a4232dc92246f54fa68d9d7b",
    )
    version(
        "7.6.2",
        sha256="1ba92942aa39b49771677cc8bf1c94a0b4350eb45bf3009318a6c2350b46a276",
    )

    variant("shared", default=True, description="Build shared libraries")
    variant("openmp", default=True, description="Build with OpenMP support")
    variant("fortran", default=False, description="Build Fortran bindings")
    variant("python", default=False, description="Build Python bindings")
    variant("memory_pool", default=True, description="Build with memory pool")
    variant("elpa", default=False, description="Use ELPA")
    variant("dlaf", default=False, when="@7.5.0:", description="Use DLA-Future")
    variant("vdwxc", default=False, description="Enable libvdwxc support")
    variant("scalapack", default=False, description="Enable scalapack support")
    variant("magma", default=False, description="Enable MAGMA support")
    variant(
        "nlcglib", default=False, description="Enable robust wave function optimization"
    )
    variant("wannier90", default=False, description="Enable Wannier90 library")
    variant(
        "build_type",
        default="Release",
        description="CMake build type",
        values=("Debug", "Release", "RelWithDebInfo"),
    )
    variant("apps", default=True, description="Build applications")
    variant("tests", default=False, description="Build tests")
    variant(
        "single_precision",
        default=False,
        description="Use single precision arithmetics",
    )
    variant(
        "profiler",
        default=True,
        description="Use internal profiler to measure execution time",
    )
    variant("nvtx", default=False, description="Use NVTX profiler")

    variant(
        "pugixml",
        default=False,
        description="Enable direct reading of UPF v2 pseudopotentials",
    )

    conflicts("+tests~pugixml")
    depends_on("pugixml", when="+pugixml")

    depends_on("cxx", type="build")
    depends_on("c", type="build")
    depends_on("fortran", type="build")

    depends_on("cmake@3.23:", type="build")
    depends_on("mpi")
    depends_on("gsl")
    depends_on("blas")
    depends_on("lapack")
    depends_on("fftw-api@3")
    depends_on("libxc@4.0.0:")
    depends_on("spglib")
    depends_on("hdf5+hl")
    depends_on("pkgconfig", type="build")

    # Python module
    depends_on("python", when="+python", type=("build", "run"))
    depends_on("python", when="@:6", type=("build", "run"))
    depends_on("py-numpy", when="+python", type=("build", "run"))
    depends_on("py-scipy", when="+python", type=("build", "run"))
    depends_on("py-h5py", when="+python", type=("build", "run"))
    depends_on("py-mpi4py", when="+python", type=("build", "run"))
    depends_on("py-pyyaml", when="+python", type=("build", "run"))
    depends_on("py-mpi4py", when="+python", type=("build", "run"))
    depends_on("py-voluptuous", when="+python", type=("build", "run"))
    depends_on("py-pybind11", when="+python", type=("build", "run"))
    extends("python", when="+python")

    depends_on("magma", when="+magma")

    depends_on("nlcglib", when="+nlcglib")
    depends_on("nlcglib+rocm", when="+nlcglib+rocm")
    depends_on("nlcglib+cuda", when="+nlcglib+cuda")

    depends_on("libvdwxc@0.3.0:+mpi", when="+vdwxc")

    depends_on("scalapack", when="+scalapack")

    with when("+dlaf"):
        depends_on("dla-future@0.3.0:")
        depends_on("dla-future +scalapack", when="+scalapack")
        depends_on("dla-future +cuda", when="+cuda")
        depends_on("dla-future +rocm", when="+rocm")

        conflicts("^pika@:0.22.1", when="+cuda")
        conflicts("^pika@:0.22.1", when="+rocm")

    depends_on("rocblas", when="+rocm")
    depends_on("rocsolver", when="@7.5.0: +rocm")

    conflicts("^libxc@5.0.0")  # known to produce incorrect results
    conflicts("+single_precision")
    conflicts("+scalapack", when="^cray-libsci")

    # Propagate openmp to blas
    depends_on(
        "openblas threads=openmp", when="+openmp ^[virtuals=blas,lapack] openblas"
    )
    depends_on("amdblis threads=openmp", when="+openmp ^[virtuals=blas] amdblis")
    depends_on("blis threads=openmp", when="+openmp ^[virtuals=blas] blis")
    depends_on(
        "intel-oneapi-mkl threads=openmp",
        when="+openmp ^[virtuals=blas,lapack,fftw-api] intel-oneapi-mkl",
    )
    depends_on(
        "intel-oneapi-mkl+cluster",
        when="+scalapack ^[virtuals=blas,lapack,fftw-api] intel-oneapi-mkl",
    )

    # MKLConfig.cmake introduced in 2021.3
    conflicts("intel-oneapi-mkl@:2021.2", when="^intel-oneapi-mkl")

    depends_on("wannier90", when="+wannier90")
    depends_on("wannier90+shared", when="+wannier90+shared")

    depends_on("elpa+openmp", when="+elpa+openmp")
    depends_on("elpa~openmp", when="+elpa~openmp")

    depends_on("eigen@3.4.0:", when="+tests")

    depends_on("costa+shared")

    with when("+memory_pool"):
        depends_on("umpire~cuda~rocm", when="~cuda~rocm")
        depends_on("umpire+cuda~device_alloc", when="+cuda")
        depends_on("umpire+rocm~device_alloc", when="+rocm")

    patch("fj.patch", when="@7.6.2: %fj")

    def cmake_args(self):
        spec = self.spec

        cm_label = "SIRIUS_"

        args = [
            self.define_from_variant(cm_label + "USE_OPENMP", "openmp"),
            self.define_from_variant(cm_label + "USE_ELPA", "elpa"),
            self.define_from_variant(cm_label + "USE_MAGMA", "magma"),
            self.define_from_variant(cm_label + "USE_NLCGLIB", "nlcglib"),
            self.define_from_variant(cm_label + "USE_VDWXC", "vdwxc"),
            self.define_from_variant(cm_label + "USE_MEMORY_POOL", "memory_pool"),
            self.define_from_variant(cm_label + "USE_SCALAPACK", "scalapack"),
            self.define_from_variant(cm_label + "USE_DLAF", "dlaf"),
            self.define_from_variant(cm_label + "CREATE_FORTRAN_BINDINGS", "fortran"),
            self.define_from_variant(cm_label + "CREATE_PYTHON_MODULE", "python"),
            self.define_from_variant(cm_label + "USE_CUDA", "cuda"),
            self.define_from_variant(cm_label + "USE_ROCM", "rocm"),
            self.define_from_variant(cm_label + "BUILD_APPS", "apps"),
            self.define_from_variant(cm_label + "USE_FP32", "single_precision"),
            self.define_from_variant(cm_label + "USE_PROFILER", "profiler"),
            self.define_from_variant(cm_label + "USE_NVTX", "nvtx"),
            self.define_from_variant(cm_label + "USE_WANNIER90", "wannier90"),
            self.define_from_variant(cm_label + "USE_PUGIXML", "pugixml"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("BUILD_TESTING", "tests"),
        ]

        lapack = spec["lapack"]
        blas = spec["blas"]

        args.extend(
            [
                self.define("LAPACK_FOUND", "true"),
                self.define("LAPACK_LIBRARIES", lapack.libs.joined(";")),
                self.define("BLAS_FOUND", "true"),
                self.define("BLAS_LIBRARIES", blas.libs.joined(";")),
            ]
        )

        if "+scalapack" in spec and "^cray-libsci" not in spec:
            args.extend(
                [
                    self.define(cm_label + "SCALAPACK_FOUND", "true"),
                    self.define(
                        cm_label + "SCALAPACK_INCLUDE_DIRS",
                        spec["scalapack"].prefix.include,
                    ),
                    self.define(
                        cm_label + "SCALAPACK_LIBRARIES",
                        spec["scalapack"].libs.joined(";"),
                    ),
                ]
            )

        if "^cray-libsci" in spec:
            args.append(self.define(cm_label + "USE_CRAY_LIBSCI", "ON"))

        if spec.satisfies("^[virtuals=blas] intel-oneapi-mkl"):
            args.append(self.define(cm_label + "USE_MKL", "ON"))

            if spec.satisfies("@7.6.0:"):
                mkl_mapper = {
                    "threading": {
                        "none": "sequential",
                        "openmp": "gnu_thread",
                        "tbb": "tbb_thread",
                    },
                    "mpi": {
                        "intel-oneapi-mpi": "intelmpi",
                        "mpich": "mpich",
                        "openmpi": "openmpi",
                    },
                }

                mkl_threads = mkl_mapper["threading"][
                    spec["intel-oneapi-mkl"].variants["threads"].value
                ]

                mpi_provider = spec["mpi"].name
                if mpi_provider in ["mpich", "cray-mpich", "mvapich", "mvapich2"]:
                    mkl_mpi = mkl_mapper["mpi"]["mpich"]
                else:
                    mkl_mpi = mkl_mapper["mpi"][mpi_provider]

                args.extend(
                    [
                        self.define("MKL_INTERFACE", "lp64"),
                        self.define("MKL_THREADING", mkl_threads),
                        self.define("MKL_MPI", mkl_mpi),
                    ]
                )

                if "+scalapack" in self.spec:
                    # options provided by `MKLConfig.cmake`
                    args.extend(
                        [
                            self.define("ENABLE_BLACS", "On"),
                            self.define("ENABLE_SCALAPACK", "On"),
                        ]
                    )

        if "+elpa" in spec:
            elpa_incdir = os.path.join(spec["elpa"].headers.directories[0], "elpa")
            args.append(self.define(cm_label + "ELPA_INCLUDE_DIR", elpa_incdir))

        if "+cuda" in spec:
            cuda_arch = spec.variants["cuda_arch"].value
            if cuda_arch[0] != "none":
                args.append(
                    self.define("CMAKE_CUDA_ARCHITECTURES", ";".join(cuda_arch))
                )

        if "+rocm" in spec:
            archs = ",".join(self.spec.variants["amdgpu_target"].value)
            args.extend([self.define("CMAKE_HIP_ARCHITECTURES", archs)])

        return args
