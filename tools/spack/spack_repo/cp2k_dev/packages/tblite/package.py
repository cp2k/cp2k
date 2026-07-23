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

    maintainers("awvwgk")

    license("LGPL-3.0-or-later")

    version("main", branch="main")
    version("0.7.0", sha256="3a7cb4602101e828caf41c38ca5e30f82de82d0d26d5db40168acdcad3462b92")
    version("0.6.0", sha256="372281aedb89234168d00eb691addb303197a9462a9c55d145c835f2cf5e8b42")
    version("0.5.0", sha256="e8a70b72ed0a0db0621c7958c63667a9cd008c97c868a4a417ff1bc262052ea8")
    version("0.4.0", sha256="5c2249b568bfd3b987d3b28f2cbfddd5c37f675b646e17c1e750428380af464b")
    version("0.3.0", sha256="46d77c120501ac55ed6a64dea8778d6593b26fb0653c591f8e8c985e35884f0a")

    variant("openmp", default=True, description="Use OpenMP parallelisation")
    variant("trexio", default=False, description="Enable TREXIO support", when="@0.7.0:")
    variant("hdf5", default=False, description="Enable HDF5 support", when="@0.7.0:")

    depends_on("c", type="build")  # generated
    depends_on("fortran", type="build")  # generated

    depends_on("blas")
    depends_on("lapack")
    depends_on("trexio", when="+trexio")
    depends_on("hdf5", when="+hdf5")


class CMakeBuilder(cmake.CMakeBuilder):
    def cmake_args(self):
        args = [self.define_from_variant("WITH_OpenMP", "openmp")]
        if self.spec.satisfies("@0.7.0:"):
            args += [
                self.define_from_variant("TBLITE_WITH_TREXIO", "trexio"),
                self.define_from_variant("TBLITE_WITH_HDF5", "hdf5"),
            ]
        return args
