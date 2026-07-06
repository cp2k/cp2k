# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.makefile import MakefilePackage

from spack.package import *


class Libgint(MakefilePackage):
    """libGint: GPU library for four-centre two-electron integrals and
    Hartree–Fock exchange calculations."""

    homepage = "https://github.com/marcelloPuligheddu/libGint"
    git = "https://github.com/marcelloPuligheddu/libGint.git"
    url = "https://github.com/marcelloPuligheddu/libGint/archive/refs/tags/release_v1.tar.gz"

    maintainers = ["marcelloPuligheddu", "mkrack"]

    license("MIT", checked_by="mkrack")

    version(
        "release_v1",
        url="https://github.com/marcelloPuligheddu/libGint/archive/refs/tags/release_v1.tar.gz",
        sha256="cc0dfeb6022ebfe0c3028a131045ff49b4c1005ad9cb77c5736fb3dac045b192",
    )

    depends_on("cuda")

    def install(self, spec, prefix):
        make("install", "PREFIX={0}".format(prefix))
        install_tree("src", prefix.include)

    @property
    def libs(self):
        return find_libraries(["libcp2kGint"], root=self.prefix, recursive=True, shared=False)
