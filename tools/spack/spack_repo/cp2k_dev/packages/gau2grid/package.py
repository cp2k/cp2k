# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator

from spack.package import *


class Gau2grid(CMakePackage):
    """Gau2Grid is a library for efficient evaluation of Gaussian
    basis functions on grids."""

    homepage = "https://gau2grid.readthedocs.io"
    git = "https://github.com/psi4/gau2grid.git"
    url = "https://github.com/psi4/gau2grid/archive/refs/tags/v2.0.8.tar.gz"

    maintainers("awvwgk", "mkrack")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("2.0.9", sha256="7879bdddf3a52cd2a051086215977822bbe8d1af927fcf5b4fb0256a38b8a76c")
    version("2.0.8", sha256="c5f445344a465c1d9afc6516544dc4a2fba588af7ba0f1ac1a6b538260f0cd96")

    depends_on("cmake@3.12:", type="build")
    depends_on("ninja@1.10:", type="build")
    depends_on("python@2.12:", type="build")
    depends_on("py-numpy", type="build")

    generator("ninja")

    def cmake_args(self):
        args = [
            self.define("MAX_AM", 8),
            self.define("ENABLE_XHOST", False),
            self.define("INSTALL_PYMOD", False),
        ]
        return args
