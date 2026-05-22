# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator

from spack.package import *


class Integratorxx(CMakePackage):
    """Reuseable density functional theory (DFT) grid library."""

    homepage = "https://github.com/wavefunction91/IntegratorXX"
    git = "https://github.com/wavefunction91/IntegratorXX.git"
    url = "https://github.com/wavefunction91/IntegratorXX/archive/refs/tags/v1.0.0.tar.gz"

    maintainers("awvwgk")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("1.0.0", sha256="d2826439d14b3f716ffd57a07d1d407de029a80c0e0446998b8f1339d5085b9c")

    depends_on("cxx", type="build")
    depends_on("cmake@3.20:", type="build")
    depends_on("ninja@1.10:", type="build")

    generator("ninja")
