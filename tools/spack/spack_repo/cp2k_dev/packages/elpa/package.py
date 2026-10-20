# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.elpa.package import Elpa as BuiltinElpa

from spack.package import *


class Elpa(BuiltinElpa):
    # GCC 16 honors the module's private default for BIND(C) symbols.
    patch("public-sirius-c-api.patch", when="@2026.02.002")
    # Backports the width-1 tail block fixes from ELPA branch bugfix_elpa2_65_65_64.
    patch("elpa-2026.02.002-nx1-tail-block.patch", when="@2026.02.002")
