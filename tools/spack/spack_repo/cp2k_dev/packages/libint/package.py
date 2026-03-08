# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libint(CMakePackage):
    """Libint is a high-performance library for computing Gaussian integrals in quantum mechanics.
    This recipe is tailored for CP2K and requires a pre-built libint archive with the source code
    generated for a specific maximum l value.
    """

    homepage = "https://github.com/evaleev/libint"
    url = "https://www.cp2k.org/static/downloads/libint-v2.13.1.tar.xz"

    maintainers("mkrack")

    license("LGPL-3.0-only")

    version(
        "2.13.1-cp2k-lmax-7",
        sha256="ba19571deefa5c3063e620210c87cf706ef8d75b47f3b0d66547f889f906b3dd",
    )
    version(
        "2.13.1-cp2k-lmax-6",
        sha256="a7990c4862d549e6328596f20b8a7ece40ded396b953d5d5d7c3200f13e37429",
    )
    version(
        "2.13.1-cp2k-lmax-5",
        sha256="527000f915bea9879391273b96689a8c2b98beeec37cb64637e3c11dd421a5e8",
    )
    version(
        "2.13.1-cp2k-lmax-4",
        sha256="8ff388fbf171635420fdfdbafc7dc949e9cf4c0b6a62f23dac3af8b1c942d407",
    )

    variant("fortran", default=True, description="Build Fortran interface")
    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=False, description="Build shared library")

    # Build dependencies
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build", when="+fortran")

    depends_on("boost", type="build")
    depends_on("eigen", type="build")

    def setup_build_environment(self, env):
        env.append_flags("CXXFLAGS", "-g1")

    def cmake_args(self):
        args = [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define_from_variant("LIBINT2_ENABLE_FORTRAN", "fortran"),
        ]
        return args

    @property
    def libs(self):
        return find_libraries("libint2", self.spec.prefix, shared=True, recursive=True)
