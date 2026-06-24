# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator

from spack.package import *


class Integratorxx(CMakePackage):
    """IntegratorXX is a modern C++ library for the generation of atomic and molecular grids for
    quantum chemistry calculations. Among the most important applications of these grids is the
    evaluation of exchange--correlation (XC) related quantities (energies, potentials, etc)
    required for density functional theory calculations.
    """

    homepage = "https://github.com/wavefunction91/IntegratorXX"
    git = "https://github.com/wavefunction91/IntegratorXX.git"
    url = "https://github.com/wavefunction91/IntegratorXX/archive/refs/tags/v1.0.0.tar.gz"

    maintainers("awvwgk", "mkrack")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("1.0.0", sha256="d2826439d14b3f716ffd57a07d1d407de029a80c0e0446998b8f1339d5085b9c")

    variant(
        "header_only", default=False, description="Force header-only build (no compiled library)"
    )
    variant("pic", default=True, description="Build position independent code")

    depends_on("cxx", type="build")
    depends_on("cmake@3.17:", type="build")
    depends_on("ninja@1.10:", type="build")

    generator("ninja")

    def cmake_args(self):
        args = [
            self.define_from_variant("INTEGRATORXX_HEADER_ONLY", "header_only"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
        ]
        return args
