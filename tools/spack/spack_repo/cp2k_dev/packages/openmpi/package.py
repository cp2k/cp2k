# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.openmpi.package import Openmpi as BuiltinOpenmpi

from spack.package import *


class Openmpi(BuiltinOpenmpi):
    # Avoid orphaned FIN/ACK packets when the shared-memory BTL goes idle.
    patch("pml-ob1-pending.patch", when="@5.0.10")
