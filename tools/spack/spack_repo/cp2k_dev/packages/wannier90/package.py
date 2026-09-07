# author: Thomas D. Kuehne, tkuehne@cp2k.org

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Wannier90(CMakePackage):
    """Wannier90 v4 with the CMake-exported Fortran library interface."""

    homepage = "https://wannier.org"
    url = "https://github.com/wannier-developers/wannier90/archive/refs/tags/v4.0.2.tar.gz"

    license("LGPL-2.1-or-later")

    version("4.0.2", sha256="2d48b371eefa8b58a6c8088c1bdffc13fe3e761111e15c8566e2ee055d8bcdb0")

    variant("mpi", default=True, description="Build with MPI support")
    variant("shared", default=True, description="Build a shared library")
    variant("pic", default=True, description="Build position independent code")

    depends_on("cmake@3.25:", type="build")
    depends_on("c", type="build")
    depends_on("fortran", type="build")
    depends_on("blas")
    depends_on("lapack")
    depends_on("mpi", when="+mpi")

    def cmake_args(self):
        args = [
            self.define_from_variant("WANNIER90_MPI", "mpi"),
            self.define_from_variant("WANNIER90_SHARED_LIBS", "shared"),
            self.define_from_variant("CMAKE_POSITION_INDEPENDENT_CODE", "pic"),
            self.define("WANNIER90_INSTALL", True),
            self.define("WANNIER90_TEST", False),
            self.define("BLA_SIZEOF_INTEGER", 4),
            self.define("BLAS_LIBRARIES", ";".join(self.spec["blas"].libs)),
            self.define("LAPACK_LIBRARIES", ";".join(self.spec["lapack"].libs)),
        ]
        if "+mpi" in self.spec:
            args.append(self.define("MPI_Fortran_COMPILER", self.spec["mpi"].mpifc))
        return args

    @run_after("install")
    def install_license(self):
        license_dir = join_path(self.prefix.share, "licenses", "wannier90")
        mkdirp(license_dir)
        install(join_path(self.stage.source_path, "LICENSE"), license_dir)
