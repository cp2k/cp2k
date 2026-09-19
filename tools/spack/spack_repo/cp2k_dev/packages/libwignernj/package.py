# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libwignernj(CMakePackage):
    """Exact evaluation of angular momentum coupling coefficients."""

    homepage = "https://github.com/susilehtola/libwignernj"
    url = "https://github.com/susilehtola/libwignernj/archive/refs/tags/v0.8.0.tar.gz"
    git = "https://github.com/susilehtola/libwignernj.git"

    maintainers("susilehtola")
    license("BSD-3-Clause", checked_by="susilehtola")

    version("main", branch="main")
    version("0.8.0", sha256="7220cea92652040d6456aba92ff151124d9c69ce8695840490c18dd25a0da80c")

    variant("fortran", default=False, description="Build the Fortran interface")
    variant("shared", default=True, description="Build shared libraries")

    depends_on("c", type="build")
    depends_on("fortran", type="build", when="+fortran")

    def cmake_args(self):
        return [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("WIGNERNJ_BUILD_FORTRAN", "fortran"),
            self.define("CMAKE_POSITION_INDEPENDENT_CODE", True),
            self.define("WIGNERNJ_BUILD_TESTS", self.run_tests),
            self.define("WIGNERNJ_BUILD_CXX_TESTS", False),
            self.define("WIGNERNJ_BUILD_EXAMPLES", False),
            self.define_from_variant("WIGNERNJ_BUILD_LTO", "ipo"),
        ]
