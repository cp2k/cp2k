# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.elpa.package import Elpa as BuiltinElpa

from spack.package import *


class Elpa(BuiltinElpa):
    # GCC 16 honors the module's private default for BIND(C) symbols.
    patch("public-sirius-c-api.patch", when="@2026.02.002")
