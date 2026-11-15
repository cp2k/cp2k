# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.libfci.package import Libfci as BuiltinLibfci

from spack.package import *


class Libfci(BuiltinLibfci):
    # Plain char is unsigned on aarch64, which breaks the fermionic signs.
    patch("libfci-0.1.0-signed-char.patch", when="@0.1.0")
