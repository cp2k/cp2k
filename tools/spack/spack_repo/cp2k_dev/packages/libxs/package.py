# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.libxs.package import Libxs as BuiltinLibxs

from spack.package import *


class Libxs(BuiltinLibxs):
    patch("libxs-1.0.0-jit-handle.patch", when="@1.0.0")
