# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems import cmake
from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Tblite(CMakePackage):
    """Light-weight tight-binding framework"""

    homepage = "https://tblite.readthedocs.io"
    url = "https://github.com/tblite/tblite/releases/download/v0.3.0/tblite-0.3.0.tar.xz"
    git = "https://github.com/tblite/tblite.git"

    maintainers("awvwgk", "mkrack")

    license("LGPL-3.0-or-later")

    version("0.6.0", sha256="372281aedb89234168d00eb691addb303197a9462a9c55d145c835f2cf5e8b42")
    version("0.5.0", sha256="e8a70b72ed0a0db0621c7958c63667a9cd008c97c868a4a417ff1bc262052ea8")
    version("0.4.0", sha256="5c2249b568bfd3b987d3b28f2cbfddd5c37f675b646e17c1e750428380af464b")
    version("0.3.0", sha256="46d77c120501ac55ed6a64dea8778d6593b26fb0653c591f8e8c985e35884f0a")

    build_system("cmake")

    variant("openmp", default=True, description="Use OpenMP parallelisation")

    depends_on("c", type="build")  # generated
    depends_on("fortran", type="build")  # generated

    depends_on("blas")
    depends_on("lapack")

    #    depends_on("mctc-lib@0.3: build_system=cmake")
    #    depends_on("simple-dftd3@0.3: build_system=cmake")
    #    depends_on("dftd4@3: build_system=cmake")
    #    depends_on("toml-f build_system=cmake")

    depends_on("pkgconfig", type="build")


class CMakeBuilder(cmake.CMakeBuilder):
    def cmake_args(self):
        return [self.define_from_variant("WITH_OpenMP", "openmp")]
