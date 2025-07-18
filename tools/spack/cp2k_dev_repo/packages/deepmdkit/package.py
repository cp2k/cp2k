# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
from spack_repo.builtin.build_systems.cmake import CMakePackage


class Deepmdkit(CMakePackage):
    """DeePMD-kit is a package written in Python/C++, designed to minimize the effort
    required to build deep learning-based model of interatomic potential energy
    and force field and to perform molecular dynamics (MD). This brings new hopes
    to addressing the accuracy-versus-efficiency dilemma in molecular simulations.
    Applications of DeePMD-kit span from finite molecules to extended systems and
    from metallic systems to chemically bonded systems."""

    homepage = "https://github.com/deepmodeling/deepmd-kit"
    url = "https://github.com/deepmodeling/deepmd-kit/archive/v3.0.2.tar.gz"
    git = "https://github.com/deepmodeling/deepmd-kit.git"

    license("LGPL-3.0", checked_by="mkrack")

    maintainers("mkrack")

    version("3.1.0", sha256="45f13df9ed011438d139a7f61416b8d7940f63c47fcde53180bfccd60c9d22ee")
    version("3.0.2", sha256="b828d3a44730ea852505abbdb24ea5b556f2bf8b16de5a9c76018ed1ced7121b")

    variant("cxx_standard", default="17", description="Required CXX standard")
    variant("cxx_standard_required", default=True, description="Require a specific CXX standard")

    variant("enable_tensorflow", default=False, description="Enable TensorFlow interface")
    variant("enable_pytorch", default=False, description="Enable PyTorch interface")
    variant("enable_jax", default=False, description="Enable JAX interface")
    variant("build_testing", default=False, description="Build test and enable coverage")
    variant("enable_native_optimization", default=False, description="Enable native optimization")
    variant("dp_using_c_api", default=True, description="Build third-party interface with C API")

    depends_on("c", type="build")
    depends_on("cxx", type="build")

    depends_on("py-torch", when="+enable_pytorch")

    root_cmakelists_dir = "source"

    build_targets = ["deepmd_c"]

    def cmake_args(self):
        args = [
            self.define_from_variant("CXX_STANDARD", "cxx_standard"),
            self.define_from_variant("CXX_STANDARD_REQUIRED", "cxx_standard_required"),
            self.define_from_variant("ENABLE_TENSORFLOW", "enable_tensorflow"),
            self.define_from_variant("ENABLE_PYTORCH", "enable_pytorch"),
            self.define_from_variant("ENABLE_JAX", "enable_jax"),
            self.define_from_variant("BUILD_TESTING", "build_testing"),
            self.define_from_variant("ENABLE_NATIVE_OPTIMIZATION", "enable_native_optimization"),
            self.define_from_variant("EDP_USING_C_AP", "dp_using_c_api"),
        ]
        return args
