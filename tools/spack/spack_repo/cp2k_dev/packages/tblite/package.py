# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems import cmake, meson
from spack_repo.builtin.build_systems.cmake import CMakePackage
from spack_repo.builtin.build_systems.meson import MesonPackage

from spack.package import *


class Tblite(CMakePackage, MesonPackage):
    """Light-weight tight-binding framework"""

    homepage = "https://tblite.readthedocs.io"
    url = "https://github.com/tblite/tblite/releases/download/v0.3.0/tblite-0.3.0.tar.xz"
    git = "https://github.com/tblite/tblite.git"

    maintainers("awvwgk")

    license("LGPL-3.0-or-later")

    version("0.5.0", sha256="e8a70b72ed0a0db0621c7958c63667a9cd008c97c868a4a417ff1bc262052ea8")
    version("0.4.0", sha256="5c2249b568bfd3b987d3b28f2cbfddd5c37f675b646e17c1e750428380af464b")
    version("0.3.0", sha256="46d77c120501ac55ed6a64dea8778d6593b26fb0653c591f8e8c985e35884f0a")

    build_system("cmake", "meson", default="meson")

    variant("openmp", default=True, description="Use OpenMP parallelisation")
    variant("python", default=False, description="Build Python extension module")

    depends_on("c", type="build")  # generated
    depends_on("fortran", type="build")  # generated

    depends_on("blas")
    depends_on("lapack")

    #    for build_system in ["cmake", "meson"]:
    #        depends_on(
    #                f"mctc-lib@0.5: build_system={build_system}", when=f"build_system={build_system}"
    #                )
    #        depends_on(
    #                f"simple-dftd3@1.2: build_system={build_system}", when=f"build_system={build_system}"
    #                )
    #        depends_on(f"dftd4@3: build_system={build_system}", when=f"build_system={build_system}")
    #        depends_on(f"toml-f build_system={build_system}", when=f"build_system={build_system}")

    depends_on("meson@0.57.2:", type="build")  # mesonbuild/meson#8377
    depends_on("pkgconfig", type="build")
    depends_on("py-cffi", when="+python")
    depends_on("py-numpy", when="+python")
    depends_on("python@3.6:", when="+python")

    extends("python", when="+python")


class MesonBuilder(meson.MesonBuilder):
    def meson_args(self):
        lapack = self.spec["lapack"].libs.names[0]
        if lapack == "lapack":
            lapack = "netlib"
        elif lapack.startswith("mkl"):
            lapack = "mkl"
        elif lapack != "openblas":
            lapack = "auto"

        return [
            "-Dlapack={0}".format(lapack),
            "-Dopenmp={0}".format(str("+openmp" in self.spec).lower()),
            "-Dpython={0}".format(str("+python" in self.spec).lower()),
        ]


class CMakeBuilder(cmake.CMakeBuilder):
    def cmake_args(self):
        return [self.define_from_variant("WITH_OpenMP", "openmp")]
