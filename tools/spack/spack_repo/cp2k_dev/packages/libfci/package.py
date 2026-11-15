# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.libfci.package import Libfci as BuiltinLibfci

from spack.package import *


class Libfci(BuiltinLibfci):
    version(
        "0.1.1",
        sha256="4364282dc8043817e7404541a347350c6d5089ea91a70ad993a7b73ca7924577",
        url="https://github.com/DCM-Uni-Paderborn/libfci/releases/download/v0.1.1/libfci-0.1.1.tar.gz",
    )
    depends_on("cmake@3.24:", when="@0.1.1:", type="build")
