# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Libfci(CMakePackage):
    """libfci is a compact determinant full-CI solver library developed by Gerald Knizia, 2010-2012,
    and packaged as a small C ABI for active-space embedding workflows. The original standalone driver
    accepted FCIDUMP-style Hamiltonians. The library interface bypasses file I/O and accepts dense
    active-space one- and two-electron matrices using CP2K's FCIDUMP conventions. The exported C entry point
    is libfci_solve. A compatibility wrapper named cp2k_fci_solve is also exported for older CP2K integration
    experiments.
    """

    maintainers("mkrack")

    homepage = "https://github.com/DCM-Uni-Paderborn/libfci"
    url = "https://www.cp2k.org/static/downloads/libfci-0.1.0.tar.gz"
    git = "https://github.com/DCM-Uni-Paderborn/libfci.git"

    license("LicenseRef-Knizia-Permissive", checked_by="mkrack")

    version("main", branch="main")
    version("0.1.0", sha256="63e8bd632e55f4b99dbaa3e905cdcd3c5dbc17aefdca1391b855f34a792bb29a")

    variant("openmp", default=True, description="Build with OpenMP parallelisation support")
    variant("pic", default=True, description="Build position independent code")
    variant("shared", default=False, description="Build shared libraries (otherwise static)")

    depends_on("cmake@3.16:", type="build")
    depends_on("cxx", type="build")

    depends_on("blas")
    depends_on("lapack")

    def cmake_args(self):
        args = [
            self.define_from_variant("LIBFCI_USE_OPENMP", "openmp"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
        ]
        return args
