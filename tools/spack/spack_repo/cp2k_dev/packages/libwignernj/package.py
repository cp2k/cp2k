# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libwignernj(CMakePackage):
    """libwignernj: exact evaluation of the Wigner 3j, 6j and 9j symbols,
    Clebsch-Gordan coefficients, Racah W-coefficients, Fano X-coefficients
    and Gaunt coefficients in C99."""

    homepage = "https://github.com/susilehtola/libwignernj"
    git = "https://github.com/susilehtola/libwignernj.git"
    url = "https://github.com/susilehtola/libwignernj/archive/refs/tags/v0.8.0.tar.gz"

    maintainers = ["susilehtola"]

    license("BSD-3-Clause", checked_by="susilehtola")

    version("main", branch="main")
    version("0.8.0", sha256="7220cea92652040d6456aba92ff151124d9c69ce8695840490c18dd25a0da80c")

    variant("fortran", default=False, description="Build the Fortran 2003 interface")
    variant("shared", default=True, description="Build a shared library")

    def cmake_args(self):
        # CP2K binds to the C interface through ISO_C_BINDING and does not
        # need the Fortran interface library
        return [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("WIGNERNJ_BUILD_FORTRAN", "fortran"),
            self.define("WIGNERNJ_BUILD_TESTS", self.run_tests),
            self.define("WIGNERNJ_BUILD_CXX_TESTS", self.run_tests),
            self.define("WIGNERNJ_BUILD_EXAMPLES", False),
        ]
