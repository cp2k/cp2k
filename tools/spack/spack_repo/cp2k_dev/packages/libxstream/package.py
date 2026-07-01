# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libxstream(CMakePackage):
    """LIBXSTREAM is an OpenCL-based accelerator backend library providing
    streams, events, and device memory management for GPU offloading. It
    implements the DBCSR ACC interface, making it a drop-in OpenCL backend for
    CP2K/DBCSR.
    """

    homepage = "https://libxstream.readthedocs.io/en/latest"
    url = "https://github.com/hfp/libxstream/archive/refs/tags/1.0.0.tar.gz"
    git = "https://github.com/hfp/libxstream.git"

    maintainers("hfp", "mkrack")

    license("BSD-3-Clause", checked_by="mkrack")

    version("main", branch="main")
    version("1.0.0", sha256="44a2823b12eb58b5eaf97649244b93dcf921597ceabc718053cd28e5f59260e3")

    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=False, description="Build shared libraries (otherwise static)")

    depends_on("cmake@3.13:", type="build")
    depends_on("c", type="build")
    depends_on("cxx", type="build")

    depends_on("opencl")
    depends_on("opencl-headers")
    depends_on("opencl-icd-loader")

    def cmake_args(self):
        spec = self.spec
        args = [
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
        ]
        return args
