# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
import copy
import os
import os.path
import sys

import spack.platforms
import spack.util.environment
import spack.util.executable
from spack.build_environment import dso_suffix
from spack.package import *


class Cp2k(CudaPackage, CMakePackage, ROCmPackage):
    """CP2K is a quantum chemistry and solid state physics software package
    that can perform atomistic simulations of solid state, liquid, molecular,
    periodic, material, crystal, and biological systems
    """

    build_system("cmake", default="cmake")

    homepage = "https://www.cp2k.org"
    url = "https://github.com/cp2k/cp2k/releases/download/v3.0.0/cp2k-3.0.tar.bz2"
    git = "https://github.com/cp2k/cp2k.git"
    list_url = "https://github.com/cp2k/cp2k/releases"

    maintainers("dev-zero", "mtaillefumier")

    version("master", branch="master", submodules="True")

    variant("mpi", default=True, description="Enable MPI support")
    variant("openmp", default=True, description="Enable OpenMP support")
    variant(
        "smm",
        default="libxsmm",
        values=("libxsmm", "libsmm", "blas"),
        description="Library for small matrix multiplications",
    )
    variant("plumed", default=False, description="Enable PLUMED support")
    variant(
        "libint",
        default=True,
        description="Use libint, required for HFX (and possibly others)",
    )
    variant(
        "libxc", default=True, description="Support additional functionals via libxc"
    )
    variant(
        "pexsi",
        default=False,
        description="Enable the alternative PEXSI method for density matrix evaluation",
    )
    variant(
        "elpa",
        default=False,
        description="Enable optimised diagonalisation routines from ELPA",
        when="@6.1:",
    )
    variant(
        "sirius",
        default=False,
        description="Enable planewave electronic structure calculations via SIRIUS",
    )
    variant("cosma", default=False, description="Use COSMA for p?gemm")
    variant(
        "libvori",
        default=False,
        description="Enable support for Voronoi integration and BQB compression",
    )
    variant("spglib", default=False, description="Enable support for spglib")
    variant(
        "spla",
        default=False,
        description="Use SPLA off-loading functionality. Only relevant when CUDA or ROCM"
        " are enabled",
    )
    variant("pytorch", default=False, description="Enable libtorch support")
    variant("quip", default=False, description="Enable quip support")
    variant("mpi_f08", default=False, description="Use MPI F08 module")

    variant(
        "enable_regtests",
        default=False,
        description="Configure cp2k to run the regtests afterwards."
        " It build cp2k normally but put the executables in exe/cmake-build-* instead of the"
        " conventional location. This option is only relevant when regtests need to be run.",
    )

    with when("+cuda"):
        variant(
            "cuda_arch_35_k20x",
            default=False,
            description=(
                "CP2K (resp. DBCSR) has specific parameter sets for"
                " different GPU models. Enable this when building"
                " with cuda_arch=35 for a K20x instead of a K40"
            ),
        )
        variant(
            "cuda_fft",
            default=False,
            description=("Use CUDA also for FFTs in the PW part of CP2K"),
        )
        variant(
            "cuda_blas",
            default=False,
            when="@:7",  # req in CP2K v8+
            description=("Use CUBLAS for general matrix operations in DBCSR"),
        )

    HFX_LMAX_RANGE = range(4, 8)

    variant(
        "lmax",
        description="Maximum supported angular momentum (HFX and others)",
        default="5",
        values=[str(x) for x in HFX_LMAX_RANGE],
        multi=False,
    )

    depends_on("python", type="build")
    depends_on("python@3:", when="@8:", type="build")

    depends_on("blas")
    depends_on("lapack")
    depends_on("fftw-api@3")

    # Force openmp propagation on some providers of blas / fftw-api
    with when("+openmp"):
        depends_on("fftw+openmp", when="^fftw")
        depends_on("amdfftw+openmp", when="^amdfftw")
        depends_on("cray-fftw+openmp", when="^cray-fftw")
        depends_on("armpl-gcc threads=openmp", when="^armpl-gcc")
        depends_on("openblas threads=openmp", when="^openblas")
        # The Cray compiler wrappers will automatically add libsci_mp with
        # -fopenmp. Since CP2K unconditionally links blas/lapack/scalapack
        # we have to be consistent.
        depends_on("cray-libsci+openmp", when="^cray-libsci")

    with when("smm=libxsmm"):
        depends_on("libxsmm@1.17:~header-only", when="@9.1:")
        # require libxsmm-1.11+ since 1.10 can leak file descriptors in Fortran
        depends_on("libxsmm@1.11:~header-only", when="@:8.9")
        # use pkg-config (support added in libxsmm-1.10) to link to libxsmm
        depends_on("pkgconfig", type="build")
        # please set variants: smm=blas by configuring packages.yaml or install
        # cp2k with option smm=blas on aarch64
        conflicts("target=aarch64:", msg="libxsmm is not available on arm")

    with when("+libint"):
        # ... and in CP2K 7.0+ for linking to libint2
        depends_on("pkgconfig", type="build", when="@7.0:")
        # libint & libxc are always statically linked
        depends_on("libint@1.1.4:1.2", when="@3.0:6.9")
        for lmax in HFX_LMAX_RANGE:
            # libint2 can be linked dynamically again
            depends_on(
                "libint@2.6.0:+fortran tune=cp2k-lmax-{0}".format(lmax),
                when="@7.0: lmax={0}".format(lmax),
            )

    with when("+libxc"):
        depends_on("pkgconfig", when="@7.0:")
        depends_on("libxc@2.2.2:3", when="@:5")
        depends_on("libxc@4.0.3:4", when="@6.0:6.9")
        depends_on("libxc@4.0.3:4", when="@7.0:8.1")
        depends_on("libxc@5.1.3:5.1", when="@8.2:8")
        depends_on("libxc@5.1.7:5.1", when="@9:2022.2")
        depends_on("libxc@6.1:", when="@2023.1:")
        depends_on("libxc@6.2:", when="@2023.2:")

    with when("+spla"):
        depends_on("spla+cuda+fortran", when="+cuda")
        depends_on("spla+rocm+fortran", when="+rocm")

    with when("+mpi"):
        depends_on("mpi@2:")
        depends_on("mpi@3:", when="@2023.1:")
        depends_on("scalapack")
        depends_on("mpich+fortran", when="^mpich")

        conflicts("~mpi_f08", when="^mpich@4.1:")

    with when("+cosma"):
        depends_on("cosma+scalapack")
        depends_on("cosma@2.5.1:", when="@9:")
        depends_on("cosma@2.6.3:", when="@2023.2:")
        depends_on("cosma+cuda", when="+cuda")
        depends_on("cosma+rocm", when="+rocm")
        conflicts("~mpi")
        # COSMA support was introduced in 8+
        conflicts("@:7")

    with when("+elpa"):
        conflicts("~mpi", msg="elpa requires MPI")
        depends_on("elpa+openmp", when="+openmp")
        depends_on("elpa~openmp", when="~openmp")
        depends_on("elpa@2021.05:", when="@8.3:")
        depends_on("elpa@2021.11.001:", when="@9.1:")
        depends_on("elpa@2023.05.001:", when="@2023.2:")

    with when("+plumed"):
        depends_on("plumed+shared")
        depends_on("plumed+mpi", when="+mpi")
        depends_on("plumed~mpi", when="~mpi")

    # while we link statically against PEXSI, its own deps may be linked in
    # dynamically, therefore can't set this as pure build-type dependency.
    with when("+pexsi"):
        conflicts("~mpi", msg="pexsi requires MPI")
        depends_on("pexsi+fortran@0.9.0:0.9", when="@:4")
        depends_on("pexsi+fortran@0.10.0:", when="@5.0:")

    # only OpenMP should be consistently used, all other common things
    # like ELPA, SCALAPACK are independent and Spack will ensure that
    # a consistent/compatible combination is pulled into the dependency graph.
    with when("+sirius"):
        depends_on("sirius+fortran+shared")
        depends_on("sirius+cuda", when="+cuda")
        depends_on("sirius+rocm", when="+rocm")
        depends_on("sirius+openmp", when="+openmp")
        depends_on("sirius~openmp", when="~openmp")
        depends_on("sirius@7.0.0:7.0", when="@8:8.2")
        depends_on("sirius@7.2", when="@8.3:8.9")
        depends_on("sirius@7.3:", when="@9.1")
        depends_on("sirius@7.4:", when="@2023.2")
        conflicts("~mpi", msg="SIRIUS requires MPI")
        # sirius support was introduced in 7, but effectively usable starting from CP2K 9
        conflicts("@:8")

    with when("+libvori"):
        depends_on("libvori@201219:", when="@8.1")
        depends_on("libvori@210412:", when="@8.2:")
        depends_on("libvori@220621:", when="@2023.1:")
        # libvori support was introduced in 8+
        conflicts("@:7")

    # the bundled libcusmm uses numpy in the parameter prediction (v7+)
    # which is written using Python 3
    depends_on("py-numpy", when="@7:+cuda")
    depends_on("python@3.6:", when="@7:+cuda")
    depends_on("py-fypp")

    depends_on("spglib", when="+spglib")

    # Apparently cp2k@4.1 needs an "experimental" version of libwannier.a
    # which is only available contacting the developer directly. See INSTALL
    # in the stage of cp2k@4.1
    depends_on("wannier90", when="@3.0+mpi")

    with when("build_system=cmake"):
        depends_on("dbcsr@2.6:")
        depends_on("dbcsr+openmp", when="+openmp")
        depends_on("dbcsr+cuda", when="+cuda")
        depends_on("dbcsr+rocm", when="+rocm")

    # CP2K needs compiler specific compilation flags, e.g. optflags
    conflicts("%apple-clang")
    conflicts("%clang")
    conflicts("%nag")
    conflicts(
        "%aocc@:3.2",
        msg="Please use AOCC 4.0+ that better support modern Fortran features CP2K requires",
    )

    conflicts(
        "~openmp", when="@8:", msg="Building without OpenMP is not supported in CP2K 8+"
    )

    # We only support specific cuda_archs for which we have parameter files
    # for optimal kernels. Note that we don't override the cuda_archs property
    # from the parent class, since the parent class defines constraints for all
    # versions. Instead just mark all unsupported cuda archs as conflicting.

    supported_cuda_arch_list = ("35", "37", "60", "70", "80")
    supported_rocm_arch_list = (
        "gfx906",
        "gfx908",
        "gfx90a",
        "gfx90a:xnack-",
        "gfx90a:xnack+",
    )
    gpu_map = {
        "35": "K40",
        "37": "K80",
        "60": "P100",
        "70": "V100",
        "80": "A100",
        "gfx906": "Mi50",
        "gfx908": "Mi100",
        "gfx90a": "Mi250",
        "gfx90a:xnack-": "Mi250",
        "gfx90a:xnack+": "Mi250",
    }
    cuda_msg = "cp2k only supports cuda_arch {0}".format(supported_cuda_arch_list)
    rocm_msg = "cp2k only supports amdgpu_target {0}".format(supported_rocm_arch_list)

    conflicts("+cuda", when="cuda_arch=none")

    # ROCm already emits an error if +rocm amdgpu_target=none is given

    with when("+cuda"):
        for arch in CudaPackage.cuda_arch_values:
            if arch not in supported_cuda_arch_list:
                conflicts("+cuda", when="cuda_arch={0}".format(arch), msg=cuda_msg)

    with when("+rocm"):
        for arch in ROCmPackage.amdgpu_targets:
            if arch not in supported_rocm_arch_list:
                conflicts("+rocm", when="amdgpu_target={0}".format(arch), msg=rocm_msg)

    def cmake_args(self):
        spec = self.spec
        args = []

        if "+cuda" in spec:
            gpu_ver = gpu_map[spec.variants["cuda_arch"].value[0]]
            args += [
                self.define("CP2K_USE_ACCEL", "CUDA"),
                self.define("CP2K_WITH_GPU", gpu_ver),
            ]

        if "+rocm" in spec:
            gpu_ver = gpu_map[spec.variants["amdgpu_target"].value[0]]
            args += [
                self.define("CP2K_USE_ACCEL", "HIP"),
                self.define("CP2K_WITH_GPU", gpu_ver),
            ]

        args += [
            self.define_from_variant("CP2K_ENABLE_REGTESTS", "enable_regtests"),
            self.define_from_variant("CP2K_USE_ELPA", "elpa"),
            self.define_from_variant("CP2K_USE_LIBINT2", "libint"),
            self.define_from_variant("CP2K_USE_SIRIUS", "sirius"),
            self.define_from_variant("CP2K_USE_SPLA", "spla"),
            self.define_from_variant("CP2K_USE_COSMA", "cosma"),
            self.define_from_variant("CP2K_USE_LIBXC", "libxc"),
            self.define_from_variant("CP2K_USE_LIBTORCH", "pytorch"),
            self.define_from_variant("CP2K_USE_METIS", "pexsi"),
            self.define_from_variant("CP2K_USE_SUPERLU", "pexsi"),
            self.define_from_variant("CP2K_USE_PLUMED", "plumed"),
            self.define_from_variant("CP2K_USE_SPGLIB", "spglib"),
            self.define_from_variant("CP2K_USE_VORI", "libvori"),
            self.define_from_variant("CP2K_USE_SPLA", "spla"),
            self.define_from_variant("CP2K_USE_QUIP", "quip"),
            self.define_from_variant("CP2K_USE_MPI_F08", "mpi_f08"),
        ]

        # we force the use elpa openmp threading support. might need to be revisited though
        args += [
            self.define(
                "CP2K_ENABLE_ELPA_OPENMP_SUPPORT",
                ("+elpa +openmp" in spec) or ("^elpa +openmp" in spec),
            )
        ]

        if "spla" in spec and (spec.satisfies("+cuda") or spec.satisfies("+rocm")):
            args += ["-DCP2K_USE_SPLA_GEMM_OFFLOADING=ON"]

        args += ["-DCP2K_USE_FFTW3=ON"]

        if spec.satisfies("smm=libxsmm"):
            args += ["-DCP2K_USE_LIBXSMM=ON"]
        else:
            args += ["-DCP2K_USE_LIBXSMM=OFF"]

        lapack = spec["lapack"]
        blas = spec["blas"]

        if blas.name in ["intel-mkl", "intel-parallel-studio", "intel-oneapi-mkl"]:
            args += ["-DCP2K_BLAS_VENDOR=MKL"]
            if sys.platform == "darwin":
                args += [
                    self.define("CP2K_BLAS_VENDOR", "CUSTOM"),
                    self.define("CP2K_SCALAPACK_VENDOR", "GENERIC"),
                    self.define(
                        "CP2K_SCALAPACK_LINK_LIBRARIES",
                        spec["scalapack"].libs.joined(";"),
                    ),
                ]
            else:
                args += ["-DCP2K_SCALAPACK_VENDOR=MKL"]
        else:
            args.extend(
                [
                    self.define("CP2K_LAPACK_FOUND", True),
                    self.define("CP2K_LAPACK_LINK_LIBRARIES", lapack.libs.joined(";")),
                    self.define("CP2K_BLAS_FOUND", True),
                    self.define("CP2K_BLAS_LINK_LIBRARIES", blas.libs.joined(";")),
                    self.define("CP2K_SCALAPACK_FOUND", True),
                    self.define(
                        "CP2K_SCALAPACK_INCLUDE_DIRS", spec["scalapack"].prefix.include
                    ),
                    self.define("CP2K_BLAS_VENDOR", "CUSTOM"),
                    self.define("CP2K_SCALAPACK_VENDOR", "GENERIC"),
                    self.define(
                        "CP2K_SCALAPACK_LINK_LIBRARIES",
                        spec["scalapack"].libs.joined(";"),
                    ),
                ]
            )

        return args

    pass
