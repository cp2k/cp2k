# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.makefile import MakefilePackage
from spack_repo.builtin.build_systems.cuda import CudaPackage

from spack.package import *


class Libgint(MakefilePackage, CudaPackage):
    """libGint: GPU library for four-centre two-electron integrals and
    Hartree–Fock exchange calculations."""

    homepage = "https://github.com/marcelloPuligheddu/libGint"
    git = "https://github.com/marcelloPuligheddu/libGint.git"
    url = "https://github.com/marcelloPuligheddu/libGint/archive/refs/tags/release_v1.tar.gz"

    maintainers = ["marcelloPuligheddu", "mkrack"]

    license("MIT", checked_by="mkrack")

    version("main", branch="main")
    version(
        "release_v1",
        url="https://github.com/marcelloPuligheddu/libGint/archive/refs/tags/release_v1.tar.gz",
        sha256="cc0dfeb6022ebfe0c3028a131045ff49b4c1005ad9cb77c5736fb3dac045b192",
    )

    variant(
        "cuda_arch",
        values=CudaPackage.cuda_arch_values,
        default="70",
        multi=False,
        description="CUDA architecture",
    )

    conflicts("~cuda")

    depends_on("mpi")

    def build(self, spec, prefix):
        arch = spec.variants["cuda_arch"].value
        make(f"ARCH_NUM={arch}")

    def install(self, spec, prefix):
        mkdirp(prefix.lib)
        mkdirp(prefix.include)
        install("libcp2kGint.a", prefix.lib)
        install("libgint.mod", prefix.include)

    @property
    def libs(self):
        return find_libraries(["cp2kGint"], root=self.prefix, recursive=True, shared=False)
