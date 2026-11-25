# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Taco(Package):
    """TACO is an open-source binary format for storing molecular-dynamics
    trajectories with positions, velocities and forces. This package installs
    the C API library used by CP2K."""

    homepage = "https://uni-ms.pages.git.nrw/ag-salinga/taco/"
    url = "https://gitlab.git.nrw/uni-ms/ag-salinga/taco/-/archive/0.3.1/taco-0.3.1.tar.gz"
    git = "https://gitlab.git.nrw/uni-ms/ag-salinga/taco.git"

    maintainers("nilsholle")

    license("MIT", checked_by="nilsholle")

    version("0.3.1", sha256="77841d34550972698e93edeb2c9f4061792ec4e3ab31f9d0635a3deb3f613276")

    depends_on("c", type="build")
    depends_on("rust@1.87:", type="build")

    def setup_build_environment(self, env):
        env.append_flags("RUSTFLAGS", "-C relocation-model=pic")

    def install(self, spec, prefix):
        cargo = which("cargo", required=True)
        with working_dir("c_api"):
            cargo("build", "--release")

        lib_dir = prefix.lib
        mkdirp(prefix.include, lib_dir)
        install("c_api/taco_format_c.h", prefix.include)
        install("c_api/target/release/libtaco_format_c.a", lib_dir)
