# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libxsmm(CMakePackage):
    """LIBXSMM is high performance library for small dense and sparse linear
    algebra opertions incl. GEMM and elementwise primities often seen in deep
    learning applications. It also serves as reference implementation of Tensor
    Processing Primitives (TPP), a programming abstraction for efficient and
    portable deep learning and HPC workloads. With version 2.0, LIBXSMM focuses
    on providing a complete and architecture-portable set of TPPs (small dense
    and sparse matrix operations as well as element-wise, GEMM, and BRGEMM
    primitives) from which higher-level operators such as convolutions,
    fully-connected layers, normalization, and pooling are composed. LIBXSMM
    targets Intel Architecture with Intel SSE, Intel AVX, Intel AVX2, Intel
    AVX‑512 (with VNNI and Bfloat16), and Intel AMX (Advanced Matrix
    Extensions), AArch64 (NEON, SVE, and SME), and RISC‑V (RVV). Code generation
    is mainly based on Just‑In‑Time (JIT) code specialization for
    compiler-independent performance (matrix multiplications, matrix
    transpose/copy, sparse functionality, and tensor primitives). LIBXSMM is
    suitable for "build once and deploy everywhere", i.e., no special target
    flags are needed to exploit the available performance. Supported GEMM
    datatypes are: FP64, FP32, FP16, bfloat16, BF8, HF8, MXBF8, MXHF8, int16,
    int8, MXBF6, MXHF6, MXFP4, int4, int2 and int1. Additionally, various
    non-standard low precision combinations are supported."""

    homepage = "https://github.com/libxsmm/libxsmm"
    url = "https://github.com/libxsmm/libxsmm/archive/2.0.0.tar.gz"
    git = "https://github.com/libxsmm/libxsmm.git"

    maintainers("hfp", "mkrack")

    license("BSD-3-Clause")

    version("main", branch="main")
    version("2.0.0", sha256="7e532dc5520f864ce6d7f44f3fd50365e3edb23da97dbdc54fd53845d86a290b")

    variant("shared", default=False, description="With shared libraries.")
    variant("fortran", default=True, description="With Fortran support.")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build", when="+fortran")

    depends_on("python", type="build")
    depends_on("cmake@3.13:", type="build")

    requires("target=x86_64:", "target=aarch64:")

    def cmake_args(self):
        spec = self.spec
        return [
            self.define("BUILD_SHARED_LIBS", spec.satisfies("+shared")),
            self.define("LIBXSMM_FORTRAN", spec.satisfies("+fortran")),
        ]
