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
    url = "https://www.cp2k.org/static/downloads/libint-2.14.0.tar.xz"

    maintainers("mkrack")

    license("LGPL-3.0-only")

    version(
        "2.14.1-cp2k-lmax-7",
        sha256="aa78b4868050e3e8605e1cb3b96f5a81238ae4b24d83eb01e9845f95c8999eec",
    )
    version(
        "2.14.1-cp2k-lmax-6",
        sha256="72bda0d0eb92f63600c3fdf65806ed9b4d49a0a3c706a3576e9a2c0c2619f31c",
    )
    version(
        "2.14.1-cp2k-lmax-5",
        sha256="6e61df7a3bbc298ded738698cb735f01bf43a2d6ec996360c6cd4b31346c571a",
    )
    version(
        "2.14.1-cp2k-lmax-4",
        sha256="4e71faf2ba1e9b0242c8970c79a74d3b9fa593fdedaaceddb2e2dcd2eed6e314",
    )

    variant("fortran", default=True, description="Build Fortran interface")
    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=False, description="Build shared library")

    # Build dependencies
    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build", when="+fortran")

    depends_on("boost")
    depends_on("eigen")

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
