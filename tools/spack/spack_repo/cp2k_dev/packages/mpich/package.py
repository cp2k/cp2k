# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import re
import sys

from spack_repo.builtin.build_systems.autotools import AutotoolsPackage
from spack_repo.builtin.build_systems.cuda import CudaPackage
from spack_repo.builtin.build_systems.rocm import ROCmPackage

from spack.package import *


class MpichEnvironmentModifications(PackageBase):
    """Collects the environment modifications that are usually needed for the life-cycle of
    MPICH, and derivatives.
    """

    def setup_dependent_build_environment(
        self, env: EnvironmentModifications, dependent_spec: Spec
    ) -> None:
        dependent_module = dependent_spec.package.module
        for var_name, attr_name in (
            ("MPICH_CC", "spack_cc"),
            ("MPICH_CXX", "spack_cxx"),
            ("MPICH_FC", "spack_fc"),
            ("MPICH_F90", "spack_fc"),
            ("MPICH_F77", "spack_f77"),
        ):
            if hasattr(dependent_module, attr_name):
                env.set(var_name, getattr(dependent_module, attr_name))

    def setup_build_environment(self, env: EnvironmentModifications) -> None:
        env.unset("F90")
        env.unset("F90FLAGS")

    def setup_run_environment(self, env: EnvironmentModifications) -> None:
        self.setup_mpi_wrapper_variables(env)

    def setup_dependent_package(self, module, dependent_spec):
        spec = self.spec
        spec.mpicc = join_path(self.prefix.bin, "mpicc")
        spec.mpicxx = join_path(self.prefix.bin, "mpicxx")
        # Some derived packages define the "fortran" variant, most don't. Checking on the
        # presence of ~fortran makes us default to add fortran wrappers if the variant is
        # not declared.
        if spec.satisfies("~fortran"):
            return
        spec.mpifc = join_path(self.prefix.bin, "mpif90")
        spec.mpif77 = join_path(self.prefix.bin, "mpif77")

    def setup_mpi_wrapper_variables(self, env):
        # Because MPI implementations provide compilers, they have to add to
        # their run environments the code to make the compilers available.
        env.set("MPICC", join_path(self.prefix.bin, "mpicc"))
        env.set("MPICXX", join_path(self.prefix.bin, "mpicxx"))
        env.set("MPIF77", join_path(self.prefix.bin, "mpif77"))
        env.set("MPIF90", join_path(self.prefix.bin, "mpif90"))


class Mpich(MpichEnvironmentModifications, AutotoolsPackage, CudaPackage, ROCmPackage):
    """MPICH is a high performance and widely portable implementation of
    the Message Passing Interface (MPI) standard."""

    homepage = "https://www.mpich.org"
    url = "https://www.mpich.org/static/downloads/3.0.4/mpich-3.0.4.tar.gz"
    git = "https://github.com/pmodels/mpich.git"
    list_url = "https://www.mpich.org/static/downloads/"
    list_depth = 1

    maintainers("raffenet", "yfguo")
    tags = ["e4s"]
    executables = ["^mpichversion$"]

    keep_werror = "specific"

    license("mpich2")

    version("develop", submodules=True)
    version("5.0.1", sha256="8c1832a13ddacf071685069f5fadfd1f2877a29e1a628652892c65211b1f3327")
    version("5.0.0", sha256="e9350e32224283e95311f22134f36c98e3cd1c665d17fae20a6cc92ed3cffe11")
    version("4.3.2", sha256="47d774587a7156a53752218c811c852e70ac44db9c502dc3f399b4cb817e3818")

    variant("hwloc", default=True, description="Use external hwloc package")
    variant("hydra", default=True, description="Build the hydra process manager")
    variant("romio", default=True, description="Enable ROMIO MPI I/O implementation")
    variant("verbs", default=False, description="Build support for OpenFabrics verbs.")
    variant("slurm", default=False, description="Enable Slurm support")
    variant("wrapperrpath", default=True, description="Enable wrapper rpath")
    variant(
        "pmi",
        default="default",
        description="""PMI interface.""",
        values=("default", "pmi", "pmi2", "pmix", "cray"),
        multi=False,
    )
    variant(
        "device",
        default="ch4",
        description="""Abstract Device Interface (ADI)
implementation. The ch4 device is in experimental state for versions
before 3.4.""",
        values=("ch3", "ch4", "ch3:sock"),
        multi=False,
    )
    variant(
        "netmod",
        default="ofi",
        description="""Network module. Only single netmod builds are
supported, and netmod is ignored if device is ch3:sock.""",
        values=("tcp", "mxm", "ofi", "ucx"),
        multi=False,
    )
    variant(
        "pci",
        default=(sys.platform != "darwin"),
        description="Support analyzing devices on PCI bus",
    )
    variant(
        "libxml2",
        default=True,
        description="Use libxml2 for XML support instead of the custom "
        "minimalistic implementation",
    )
    variant("argobots", default=False, description="Enable Argobots support")
    variant("fortran", default=True, description="Enable Fortran support")

    variant(
        "vci",
        default=False,
        when="@4: device=ch4",
        description="Enable multiple VCI (virtual communication "
        "interface) critical sections to improve performance "
        "of applications that do heavy concurrent MPI"
        "communications. Set MPIR_CVAR_CH4_NUM_VCIS=<N> to "
        "enable multiple vcis at runtime.",
    )

    variant(
        "datatype-engine",
        default="auto",
        description="controls the datatype engine to use",
        values=("dataloop", "yaksa", "auto"),
        when="@3.4:",
        multi=False,
    )
    for _yaksa_cond in (
        "@4.0:4.3 device=ch4 datatype-engine=auto",
        "@4.0:4.3 device=ch4 datatype-engine=yaksa",
    ):
        with when(_yaksa_cond):
            depends_on("yaksa")
            depends_on("yaksa+cuda", when="+cuda")
            depends_on("yaksa+rocm", when="+rocm")

    variant(
        "hcoll",
        default=False,
        description="Enable support for Mellanox HCOLL accelerated collective operations library",
        when="@3.3: device=ch4 netmod=ucx",
    )

    variant("xpmem", default=False, when="@3.4:", description="Enable XPMEM support")
    variant("level_zero", default=False, description="Enable level zero support")

    conflicts("datatype-engine=yaksa", when="device=ch3")
    conflicts("datatype-engine=yaksa", when="device=ch3:sock")
    conflicts("datatype-engine=dataloop", when="+cuda")
    conflicts("datatype-engine=dataloop", when="+rocm")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type="build", when="+fortran")

    depends_on("hcoll", when="+hcoll")
    depends_on("xpmem", when="+xpmem")

    # Todo: cuda can be a conditional variant, but it does not seem to work when
    # overriding the variant from CudaPackage.
    conflicts("+cuda", when="@:3.3")
    conflicts("+cuda", when="device=ch3")
    conflicts("+cuda", when="device=ch3:sock")
    conflicts("+rocm", when="@:4.0")
    conflicts("+rocm", when="device=ch3")
    conflicts("+rocm", when="device=ch3:sock")
    conflicts("+cuda", when="+rocm", msg="CUDA must be disabled to support ROCm")

    provides("mpi@:5.0", when="@5:")
    provides("mpi@:4.0", when="@:4.3")
    provides("mpi@:3.1", when="@:3.2")
    provides("mpi@:3.0", when="@:3.1")
    provides("mpi@:2.2", when="@:1.2")
    provides("mpi@:2.1", when="@:1.1")
    provides("mpi@:2.0", when="@:1.0")

    filter_compiler_wrappers("mpicc", "mpicxx", "mpif77", "mpif90", "mpifort", relative_root="bin")

    depends_on("findutils", type="build")
    depends_on("pkgconfig", type="build")

    depends_on("hwloc@2.0.0:", when="@3.3: +hwloc")
    depends_on("hwloc@2.0.0: +cuda", when="@3.3: +cuda+hwloc")

    depends_on("libfabric", when="netmod=ofi")
    depends_on("libfabric+cuda", when="+cuda netmod=ofi")
    # The ch3 ofi netmod results in crashes with libfabric 1.7
    # See https://github.com/pmodels/mpich/issues/3665
    depends_on("libfabric@:1.6", when="device=ch3 netmod=ofi")
    depends_on("libfabric@1.5:", when="@3.4: device=ch4 netmod=ofi")

    depends_on("ucx", when="netmod=ucx")
    depends_on("ucx+cuda", when="+cuda netmod=ucx")
    depends_on("mxm", when="netmod=mxm")

    # The dependencies on libpciaccess and libxml2 come from the embedded
    # hwloc, which, before version 3.3, was used only for Hydra.
    depends_on("libpciaccess", when="@:3.2+hydra+pci")
    depends_on("libxml2", when="@:3.2+hydra+libxml2")

    # Starting with version 3.3, MPICH uses hwloc directly.
    depends_on("libpciaccess", when="@3.3:+pci")
    depends_on("libxml2", when="@3.3:+libxml2")

    # Starting with version 3.3, Hydra can use libslurm for nodelist parsing
    depends_on("slurm", when="+slurm")

    depends_on("pmix", when="pmi=pmix")

    # +argobots variant requires Argobots
    depends_on("argobots", when="+argobots")

    # building from git requires regenerating autotools files
    depends_on("automake@1.15:", when="@develop", type="build")
    depends_on("libtool@2.4.4:", when="@develop", type="build")
    depends_on("m4", when="@develop", type="build")
    depends_on("autoconf@2.67:", when="@develop", type="build")

    # building with "+hwloc' also requires regenerating autotools files
    depends_on("automake@1.15:", when="@3.3 +hwloc", type="build")
    depends_on("libtool@2.4.4:", when="@3.3 +hwloc", type="build")
    depends_on("m4", when="@3.3 +hwloc", type="build")
    depends_on("autoconf@2.67:", when="@3.3 +hwloc", type="build")

    # MPICH's Yaksa submodule requires python to configure
    depends_on("python@3.0:", when="@5:", type="build")

    depends_on("cray-pmi", when="pmi=cray")
    depends_on("oneapi-level-zero", when="+level_zero")

    conflicts("device=ch4", when="@:3.2")
    conflicts("netmod=ofi", when="@:3.1.4")
    conflicts("netmod=ucx", when="device=ch3")
    conflicts("netmod=mxm", when="device=ch4")
    conflicts("netmod=mxm", when="@:3.1.3")
    conflicts("netmod=tcp", when="device=ch4")
    conflicts("pmi=pmi2", when="device=ch3 netmod=ofi")
    conflicts("pmi=pmix", when="device=ch3")
    conflicts("pmi=pmix", when="device=ch3:sock")
    conflicts("pmi=pmix", when="+hydra")
    conflicts("pmi=cray", when="+hydra")

    # MPICH does not require libxml2 and libpciaccess for versions before 3.3
    # when ~hydra is set: prevent users from setting +libxml2 and +pci in this
    # case to avoid generating an identical MPICH installation.
    conflicts("+pci", when="@:3.2~hydra")
    conflicts("+libxml2", when="@:3.2~hydra")

    # see https://github.com/pmodels/mpich/pull/5031
    conflicts("%clang@:7", when="@3.4:3.4.1")

    @classmethod
    def determine_version(cls, exe):
        output = Executable(exe)(output=str, error=str)
        match = re.search(r"MPICH Version:\s+(\S+)", output)
        return match.group(1) if match else None

    @classmethod
    def determine_variants(cls, exes, version):
        def get_spack_compiler_spec(compiler):
            spack_compilers = find_compilers([os.path.dirname(compiler)])
            actual_compiler = None
            # check if the compiler actually matches the one we want
            for spack_compiler in spack_compilers:
                if spack_compiler.cc and spack_compiler.cc == compiler:
                    actual_compiler = spack_compiler
                    break
            return actual_compiler.spec if actual_compiler else None

        def is_enabled(text):
            if text in set(["t", "true", "enabled", "enable", "with", "yes", "1"]):
                return True
            return False

        def is_disabled(text):
            if text in set(["f", "false", "disabled", "disable", "without", "no", "0"]):
                return True
            return False

        results = []
        for exe in exes:
            variants = []
            output = Executable(exe)(output=str, error=str)
            if re.search(r"--with-hwloc(-prefix)*=embedded", output):
                variants.append("~hwloc")

            if re.search(r"--with-pm=hydra", output):
                variants.append("+hydra")
            else:
                variants.append("~hydra")

            match = re.search(r"--(\S+)-romio", output)
            if match and is_enabled(match.group(1)):
                variants.append("+romio")
            elif match and is_disabled(match.group(1)):
                variants.append("~romio")

            if re.search(r"--with-ibverbs", output):
                variants.append("+verbs")
            elif re.search(r"--without-ibverbs", output):
                variants.append("~verbs")

            match = re.search(r"--enable-wrapper-rpath=(\S+)", output)
            if match and is_enabled(match.group(1)):
                variants.append("+wrapperrpath")
            match = re.search(r"--enable-wrapper-rpath=(\S+)", output)
            if match and is_disabled(match.group(1)):
                variants.append("~wrapperrpath")

            if re.search(r"--disable-fortran", output):
                variants.append("~fortran")

            match = re.search(r"--with-slurm=(\S+)", output)
            if match and is_enabled(match.group(1)):
                variants.append("+slurm")

            if re.search(r"--enable-libxml2", output):
                variants.append("+libxml2")
            elif re.search(r"--disable-libxml2", output):
                variants.append("~libxml2")

            if re.search(r"--with-thread-package=argobots", output):
                variants.append("+argobots")

            if re.search(r"--with-pmi=default", output):
                variants.append("pmi=default")
            elif re.search(r"--with-pmi=simple", output):
                variants.append("pmi=pmi")
            elif re.search(r"--with-pmi=pmi2/simple", output):
                variants.append("pmi=pmi2")
            elif re.search(r"--with-pmi=pmi", output):
                variants.append("pmi=pmi")
            elif re.search(r"--with-pmi=pmi2", output):
                variants.append("pmi=pmi2")
            elif re.search(r"--with-pmix", output):
                variants.append("pmi=pmix")

            match = re.search(r"MPICH Device:\s+(ch3|ch4)", output)
            if match:
                variants.append("device=" + match.group(1))

            match = re.search(r"--with-device=ch.\S+(ucx|ofi|mxm|tcp)", output)
            if match:
                variants.append("netmod=" + match.group(1))

            if re.search(r"--with-hcoll", output):
                variants += "+hcoll"

            match = re.search(r"MPICH CC:\s+(\S+)", output)
            if match:
                compiler = match.group(1)
                compiler_spec = get_spack_compiler_spec(compiler)
                if compiler_spec:
                    variants.append("%" + str(compiler_spec))
            results.append(" ".join(variants))
        return results

    def flag_handler(self, name, flags):
        if name == "fflags":
            # https://bugzilla.redhat.com/show_bug.cgi?id=1795817
            # https://github.com/spack/spack/issues/17934
            # TODO: we should add the flag depending on the real Fortran compiler spec and not the
            #  toolchain spec, which might be mixed.
            if any(self.spec.satisfies(s) for s in ["%gcc@10:", "%apple-clang@11:", "%clang@11:"]):
                # Note that the flag is not needed to build the package starting version 4.1
                # (see https://github.com/pmodels/mpich/pull/5840) but we keep adding the flag here
                # to avoid its presence in the MPI compiler wrappers.
                flags.append("-fallow-argument-mismatch")

        return flags, None, None

    def setup_build_environment(self, env: EnvironmentModifications) -> None:
        MpichEnvironmentModifications.setup_build_environment(self, env)
        if "pmi=cray" in self.spec:
            env.set("CRAY_PMI_INCLUDE_OPTS", "-I" + self.spec["cray-pmi"].headers.directories[0])
            env.set("CRAY_PMI_POST_LINK_OPTS", "-L" + self.spec["cray-pmi"].libs.directories[0])

    def autoreconf(self, spec, prefix):
        """Not needed usually, configure should be already there"""
        # If configure exists nothing needs to be done
        if os.path.exists(self.configure_abs_path) and not spec.satisfies("@3.3 +hwloc"):
            return
        # Else bootstrap with autotools
        bash = which("bash", required=True)
        bash("./autogen.sh")

    def configure_args(self):
        spec = self.spec
        config_args = [
            "--disable-maintainer-mode",
            "--disable-silent-rules",
            "--enable-shared",
            "--with-pm={0}".format("hydra" if "+hydra" in spec else "no"),
            "--{0}-romio".format("enable" if "+romio" in spec else "disable"),
            "--{0}-ibverbs".format("with" if "+verbs" in spec else "without"),
            "--enable-wrapper-rpath={0}".format("no" if "~wrapperrpath" in spec else "yes"),
            "--with-yaksa={0}".format(spec["yaksa"].prefix if "^yaksa" in spec else "embedded"),
            *self.with_or_without("ze", variant="level_zero"),
        ]

        # https://github.com/pmodels/mpich/commit/bbfc4cab6ade0b75ef3803a83af1cad4a262a564
        if self.spec.satisfies("@:4.2 ~hwloc"):
            config_args += self.enable_or_disable("levelzero", variant="level_zero")

        # see https://github.com/pmodels/mpich/issues/5530
        if spec.platform == "darwin":
            config_args.append("--enable-two-level-namespace")

        # hwloc configure option changed in 4.0
        if spec.satisfies("@4.0:"):
            config_args.append(
                "--with-hwloc={0}".format(spec["hwloc"].prefix if "^hwloc" in spec else "embedded")
            )
        else:
            config_args.append(
                "--with-hwloc-prefix={0}".format(
                    spec["hwloc"].prefix if "^hwloc" in spec else "embedded"
                )
            )

        config_args.extend(self.enable_or_disable("fortran"))

        if "+slurm" in spec:
            config_args.append("--with-slurm=yes")
            config_args.append("--with-slurm-include={0}".format(spec["slurm"].prefix.include))
            config_args.append("--with-slurm-lib={0}".format(spec["slurm"].prefix.lib))
        else:
            config_args.append("--with-slurm=no")

        # PMI options changed in 4.2.0
        if spec.satisfies("@4.2:"):
            # default (no option) is to build both PMIv1 and PMI2 client interfaces
            if "pmi=pmi" in spec:
                # make PMI1 the default client interface
                config_args.append("--with-pmi=pmi")
            elif "pmi=pmi2" in spec:
                # make PMI2 the default client interface
                config_args.append("--with-pmi=pmi2")
            elif "pmi=pmix" in spec:
                # use the PMIx client interface with an external PMIx library
                config_args.append("--with-pmi=pmix")
                config_args.append(f"--with-pmix={spec['pmix'].prefix}")
            elif "pmi=cray" in spec:
                # use PMI2 interface of the Cray PMI library
                config_args.append("--with-pmi=pmi2")
                config_args.append(f"--with-pmi2={spec['cray-pmi'].prefix}")
        else:
            if "pmi=pmi" in spec:
                config_args.append("--with-pmi=simple")
            elif "pmi=pmi2" in spec:
                config_args.append("--with-pmi=pmi2/simple")
            elif "pmi=pmix" in spec:
                config_args.append(f"--with-pmix={spec['pmix'].prefix}")
            elif "pmi=cray" in spec:
                config_args.append("--with-pmi=cray")

        if "+cuda" in spec:
            config_args.append("--with-cuda={0}".format(spec["cuda"].prefix))
        elif not spec.satisfies("@3.4:3.4.3"):
            # Versions from 3.4 to 3.4.3 cannot handle --without-cuda
            # (see https://github.com/pmodels/mpich/pull/5060):
            config_args.append("--without-cuda")

        if "+rocm" in spec:
            config_args.append("--with-hip={0}".format(spec["hip"].prefix))
        else:
            config_args.append("--without-hip")

        # setup device configuration
        device_config = ""
        if "device=ch4" in spec:
            device_config = "--with-device=ch4:"
        elif "device=ch3" in spec:
            device_config = "--with-device=ch3:nemesis:"

        # Do not apply any netmod if device is ch3:sock
        if "device=ch3:sock" in spec:
            device_config = "--with-device=ch3:sock"
        elif "netmod=ucx" in spec:
            device_config += "ucx"
        elif "netmod=ofi" in spec:
            device_config += "ofi"
        elif "netmod=mxm" in spec:
            device_config += "mxm"
        elif "netmod=tcp" in spec:
            device_config += "tcp"

        config_args.append(device_config)

        # Specify libfabric or ucx path explicitly, otherwise
        # configure might fall back to an embedded version.
        if "netmod=ofi" in spec:
            config_args.append("--with-libfabric={0}".format(spec["libfabric"].prefix))
        if "netmod=ucx" in spec:
            config_args.append("--with-ucx={0}".format(spec["ucx"].prefix))

        # In other cases the argument is redundant.
        if "@:3.2+hydra" in spec or "@3.3:" in spec:
            # The root configure script passes the argument to the configure
            # scripts of all instances of hwloc (there are three copies of it:
            # for hydra, for hydra2, and for MPICH itself).
            config_args += self.enable_or_disable("libxml2")

        # If +argobots specified, add argobots option
        if "+argobots" in spec:
            config_args.append("--with-thread-package=argobots")
            config_args.append("--with-argobots=" + spec["argobots"].prefix)

        if "+vci" in spec:
            config_args.append("--enable-thread-cs=per-vci")

        if "datatype-engine=yaksa" in spec:
            config_args.append("--with-datatype-engine=yaksa")
        elif "datatype-engine=dataloop" in spec:
            config_args.append("--with-datatype-engine=dataloop")
        elif "datatype-engine=auto" in spec:
            config_args.append("--with-datatype-engine=auto")

        if "+hcoll" in spec:
            config_args.append("--with-hcoll=" + spec["hcoll"].prefix)

        if "+xpmem" in spec:
            config_args.append("--with-xpmem=" + spec["xpmem"].prefix)

        return config_args

    @run_after("install")
    def cache_test_sources(self):
        """Copy the example source files after the package is installed to an
        install test subdirectory for use during `spack test run`."""
        cache_extra_test_sources(self, ["examples", join_path("test", "mpi")])

    def mpi_launcher(self):
        """Determine the appropriate launcher."""
        commands = [
            join_path(self.spec.prefix.bin, "mpirun"),
            join_path(self.spec.prefix.bin, "mpiexec"),
        ]
        if "+slurm" in self.spec:
            commands.insert(0, join_path(self.spec["slurm"].prefix.bin))
        return which(*commands, required=True)

    def run_mpich_test(self, subdir, exe, num_procs=1):
        """Compile and run the test program."""
        path = self.test_suite.current_test_cache_dir.join(subdir)
        with working_dir(path):
            src = f"{exe}.c"
            if not os.path.isfile(src):
                raise SkipTest(f"{src} is missing")

            mpicc = which(os.environ["MPICC"], required=True)
            mpicc("-Wall", "-g", "-o", exe, src)
            if num_procs > 1:
                launcher = self.mpi_launcher()
                if launcher is not None:
                    launcher("-n", str(num_procs), exe)
                    return

            test_exe = which(exe, required=True)
            test_exe()

    def test_cpi(self):
        """build and run cpi"""
        self.run_mpich_test("examples", "cpi")

    def test_finalized(self):
        """build and run finalized"""
        self.run_mpich_test(join_path("test", "mpi", "init"), "finalized")

    def test_manyrma(self):
        """build and run manyrma"""
        self.run_mpich_test(join_path("test", "mpi", "perf"), "manyrma", 2)

    def test_sendrecv(self):
        """build and run sendrecv"""
        self.run_mpich_test(join_path("test", "mpi", "basic"), "sendrecv", 2)
