# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.packages.py_torch.package import PyTorch as BuiltinPyTorch

from spack.package import *


class PyTorch(BuiltinPyTorch):
    # PyTorch's cmake/Dependencies.cmake hardcodes the source of the fp16 target
    # to "/usr/include/fp16.h" when USE_SYSTEM_FP16=ON. Spack installs fp16 into
    # its own prefix instead, so CMake fails with "Cannot find source file".
    # Build against the vendored third_party/FP16 submodule (checked out by
    # Spack) instead.
    def setup_build_environment(self, env):
        super().setup_build_environment(env)
        env.set("USE_SYSTEM_FP16", "OFF")
