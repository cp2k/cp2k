# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack_repo.builtin.build_systems.cmake import CMakePackage, generator

from spack.package import *


class Gau2grid(CMakePackage):
    """Gau2Grid is a library for the fast computation of a Gaussian and its derivative on a grid."""

    homepage = "https://github.com/psi4/gau2grid"
    git = "https://github.com/psi4/gau2grid.git"
    url = "https://github.com/psi4/gau2grid/archive/v2.0.9.tar.gz"

    maintainers("awvwgk", "mkrack")

    license("BSD-3-Clause")

    version("master", branch="master")
    version("2.0.9", sha256="7879bdddf3a52cd2a051086215977822bbe8d1af927fcf5b4fb0256a38b8a76c")
    version("2.0.8", sha256="c5f445344a465c1d9afc6516544dc4a2fba588af7ba0f1ac1a6b538260f0cd96")

    variant("max_am", default="8", description="Maximum Gaussian angular momentum")
    variant(
        "xhost",
        default=True,
        description="Enable processor-specific optimizations (-march=native)",
    )
    variant("pic", default=True, description="Build with position independent code")
    variant("shared", default=True, description="Build shared libraries")
    variant("generic", default=False, description="Enable generic/static system linking")
    variant("no_pragma", default=False, description="Disable pragma optimizations")

    conflicts("+shared~pic", msg="Shared library requires position independent code")

    depends_on("c", type="build")
    depends_on("cmake@3.12:", type="build")
    depends_on("ninja@1.10:", type="build")
    depends_on("python@3.6:", type=("build", "link"))
    depends_on("py-numpy", type=("build", "link"))

    generator("ninja")

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define("MAX_AM", spec.variants["max_am"].value),
            self.define("ENABLE_XHOST", spec.satisfies("+xhost")),
            self.define("BUILD_FPIC", spec.satisfies("+pic")),
            self.define("BUILD_SHARED_LIBS", spec.satisfies("+shared")),
            self.define("ENABLE_GENERIC", spec.satisfies("+generic")),
            self.define("DISABLE_PRAGMA", spec.satisfies("+no_pragma")),
        ]
        return args
