# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.dftd4.package import CMakeBuilder as BuiltinCMakeBuilder
from spack_repo.builtin.packages.dftd4.package import Dftd4 as BuiltinDftd4
from spack_repo.builtin.packages.dftd4.package import MesonBuilder as BuiltinMesonBuilder

from spack.package import *


class Dftd4(BuiltinDftd4):
    depends_on("mstore build_system=cmake", when="build_system=cmake")

    # GCC 16 honors the module's private default for BIND(C) symbols.
    patch("public-numerical-hessian.patch", when="@3.5:")


class MesonBuilder(BuiltinMesonBuilder):
    pass


class CMakeBuilder(BuiltinCMakeBuilder):
    def cmake_args(self):
        return [
            self.define_from_variant("WITH_OpenMP", "openmp"),
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
        ]
