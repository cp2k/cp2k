# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.simple_dftd3.package import CMakeBuilder as BuiltinCMakeBuilder
from spack_repo.builtin.packages.simple_dftd3.package import MesonBuilder as BuiltinMesonBuilder
from spack_repo.builtin.packages.simple_dftd3.package import SimpleDftd3 as BuiltinSimpleDftd3

from spack.package import *


class SimpleDftd3(BuiltinSimpleDftd3):
    # GCC 16 honors the module's private default for BIND(C) symbols.
    patch("public-c-api.patch", when="@1.4.0")


class MesonBuilder(BuiltinMesonBuilder):
    pass


class CMakeBuilder(BuiltinCMakeBuilder):
    pass
