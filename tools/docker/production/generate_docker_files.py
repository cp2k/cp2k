#!/usr/bin/env python3

# Author: Matthias Krack (October 16, 2023)

from pathlib import Path
from typing import Any
import argparse
import io
import os

# ------------------------------------------------------------------------------

cp2k_release_list = ["master", "2023.2"]  # append new releases to list
mpi_implementation_list = ["intelmpi", "mpich", "openmpi"]
target_cpu_list = ["generic", "haswell", "native", "skylake-avx512", "znver2", "znver3"]
target_gpu_list = ["P100", "V100", "A100"]  # append new GPUs to list

# ------------------------------------------------------------------------------


def main() -> None:
    mpi_choices = ["all"] + mpi_implementation_list
    release_choices = ["all"] + cp2k_release_list
    target_cpu_choices = ["all"] + target_cpu_list
    target_gpu_choices = ["all"] + target_gpu_list
    version_choices = ["psmp", "ssmp", "pdbg", "sdbg"]

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        dest="check",
        help="Check consistency with generator script",
    )
    parser.add_argument(
        "--mpi",
        choices=mpi_choices,
        default=mpi_choices[0],
        dest="mpi_implementation",
        help=(
            "Select a MPI implementation (default is to generate docker "
            f"containers for {mpi_choices[0]})"
        ),
        type=str,
    )
    parser.add_argument(
        "-j",
        "--ncores",
        default=8,
        dest="ncores",
        help=(
            "Select the number of CPU cores used for building the container "
            "and running the regression tests (default is 8)"
        ),
        type=check_ncores,
    )
    parser.add_argument(
        "--no-tests",
        action="store_true",
        dest="no_tests",
        help="The container will not include the files for regression testing",
    )
    parser.add_argument(
        "--release",
        choices=release_choices,
        default=release_choices[0],
        dest="release",
        help=(
            "Specify the CP2K release for which the docker files are generated "
            f"(default is {release_choices[0]})"
        ),
        type=str,
    )
    parser.add_argument(
        "--target-cpu",
        choices=target_cpu_choices,
        default=target_cpu_choices[0],
        dest="target_cpu",
        help=(
            "Specify the target CPU for which the docker files are generated "
            f"(default is {target_cpu_choices[0]})"
        ),
        type=str,
    )
    parser.add_argument(
        "--target-gpu",
        choices=target_gpu_choices,
        default=target_gpu_choices[0],
        dest="target_gpu",
        help=(
            "Specify the target GPU for which the docker files are generated "
            f"(default is {target_gpu_choices[0]})"
        ),
        type=str,
    )
    parser.add_argument(
        "--test",
        "--test-build",
        action="store_true",
        dest="test_build",
        help="Run a full regression test during the build step",
    )
    parser.add_argument(
        "--version",
        choices=version_choices,
        default=version_choices[0],
        dest="version",
        help=(
            "Specify the version type of the CP2K binary "
            f"(default is {version_choices[0]})"
        ),
        type=str,
    )
    args = parser.parse_args()

    version = args.version
    ncores = args.ncores
    no_tests = args.no_tests
    omp_stacksize = "16M"
    test_build = args.test_build

    if ncores > os.cpu_count():
        print(
            "WARNING: More CPU cores requested for build than available "
            f"({ncores} > {os.cpu_count()})"
        )

    for release in cp2k_release_list:
        if args.release != "all" and args.release != release:
            continue

        for mpi_implementation in mpi_implementation_list:
            if (
                args.mpi_implementation != "all"
                and args.mpi_implementation != mpi_implementation
            ):
                continue

            if mpi_implementation == "intelmpi":
                base_system = "intel/oneapi-hpckit:2023.2.1-devel-ubuntu22.04"
                if release != "master":
                    if float(release) <= 2023.2:
                        continue
            else:
                base_system = "ubuntu:22.04"

            for target_cpu in target_cpu_list:
                if args.target_cpu != "all" and args.target_cpu != target_cpu:
                    continue
                if "znver" in target_cpu and mpi_implementation == "intelmpi":
                    continue
                if release != "master":
                    if float(release) <= 2023.2 and "znver" in target_cpu:
                        continue
                name = f"{release}_{mpi_implementation}_{target_cpu}_{version}"
                with OutputFile(f"Dockerfile.{name}", args.check) as f:
                    f.write(
                        write_definition_file(
                            base_system=base_system,
                            name=name,
                            release=release,
                            version=version,
                            mpi_implementation=mpi_implementation,
                            ncores=ncores,
                            no_tests=no_tests,
                            omp_stacksize=omp_stacksize,
                            target_cpu=target_cpu,
                            target_gpu="",
                            test_build=test_build,
                        )
                    )

    # Generate docker files for CUDA
    base_system = "nvidia/cuda:12.2.0-devel-ubuntu22.04"
    for release in cp2k_release_list:
        if args.release != "all" and args.release != release:
            continue
        for mpi_implementation in mpi_implementation_list:
            if (
                args.mpi_implementation != "all"
                and args.mpi_implementation != mpi_implementation
            ) or mpi_implementation == "intelmpi":
                continue
            for target_cpu in target_cpu_list:
                if args.target_cpu != "all" and args.target_cpu != target_cpu:
                    continue
                # Restrict docker file generation for CUDA
                if target_cpu not in ["generic", "native", args.target_cpu]:
                    continue
                for igpu, target_gpu in enumerate(target_gpu_list):
                    if args.target_gpu != "all" and args.target_gpu != target_gpu:
                        continue
                    if release != "master":
                        if float(release) <= 2023.2 and igpu > 1:
                            continue
                    name = (
                        f"{release}_{mpi_implementation}_{target_cpu}"
                        f"_cuda_{target_gpu}_{version}"
                    )
                    with OutputFile(f"Dockerfile.{name}", args.check) as f:
                        f.write(
                            write_definition_file(
                                base_system=base_system,
                                name=name,
                                release=release,
                                version=version,
                                mpi_implementation=mpi_implementation,
                                ncores=ncores,
                                no_tests=no_tests,
                                omp_stacksize=omp_stacksize,
                                target_cpu=target_cpu,
                                target_gpu=target_gpu,
                                test_build=test_build,
                            )
                        )


# ------------------------------------------------------------------------------


def check_ncores(value: str) -> int:
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"{value} is an invalid number of CPU cores")
    return ivalue


# ------------------------------------------------------------------------------


def write_definition_file(
    base_system: str,
    name: str,
    release: str,
    version: str,
    mpi_implementation: str,
    ncores: int,
    no_tests: bool,
    omp_stacksize: str,
    target_cpu: str,
    target_gpu: str,
    test_build: bool,
) -> str:
    do_regtest = "/opt/cp2k/tests/do_regtest.py"

    if release == "master":
        branch = ""
        tagname = name.replace("master", r"master$(date +%Y%m%d)")
    else:
        branch = f" -b support/v{release}"
        tagname = name
        # The location of the regression test script has changed only after v2023.2
        if float(release) <= 2023.2:
            do_regtest = "/opt/cp2k/tools/regtesting/do_regtest.py"

    if test_build:
        make_target = " test"
    else:
        make_target = ""

    # Required packages for the final container
    required_packages = "g++ gcc gfortran"
    if target_gpu:
        cuda_path = "/usr/local/cuda"
        additional_exports = f"""\
export CUDA_CACHE_DISABLE=1\\n\\
export CUDA_PATH={cuda_path}\\n\\
export LD_LIBRARY_PATH=\${{LD_LIBRARY_PATH}}:\${{CUDA_PATH}}/lib64\\n\\\
"""
        arch = "local_cuda"
        cuda_packages = " libtool libtool-bin"
        cuda_environment = f"""\
# Setup CUDA environment
ENV CUDA_PATH {cuda_path}
ENV LD_LIBRARY_PATH {cuda_path}/lib64

# Disable JIT cache as there seems to be an issue with file locking on overlayfs
# See also https://github.com/cp2k/cp2k/pull/2337
ENV CUDA_CACHE_DISABLE 1
"""
        with_gpu_line = f"--enable-cuda=yes --gpu-ver={target_gpu} --with-libtorch=no"
    else:
        additional_exports = "\\"
        arch = "local"
        cuda_packages = ""
        cuda_environment = ""
        with_gpu_line = "--enable-cuda=no"

    # Default options for the regression tests
    testopts = f"--maxtasks {ncores}" " --workbasedir /mnt"

    permissive = ""
    with_mpi_line = f"--with-{mpi_implementation}=system"
    if mpi_implementation == "mpich":
        required_packages += " libmpich-dev mpich openssh-client"
        with_gcc_line = "--with-gcc=system"
    elif mpi_implementation == "openmpi":
        additional_exports += """
export OMPI_ALLOW_RUN_AS_ROOT=1\\n\\
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\\n\\
export OMPI_MCA_btl_vader_single_copy_mechanism=none\\n\\\
"""
        # SuperLU installation fails currently with system OpenMPI
        # required_packages += " libopenmpi-dev openmpi-bin openssh-client"
        with_mpi_line = f"--with-{mpi_implementation}=install"
        required_packages += " openssh-client"
        testopts = '--mpiexec \\"mpiexec --bind-to none\\" ' + testopts
        with_gcc_line = "--with-gcc=system"
    elif mpi_implementation == "intelmpi":
        permissive = "; true"
        with_gcc_line = (
            "--with-intel=system --with-intelmpi=system "
            "--with-libtorch=no --with-mkl=system"
        )

    if no_tests:
        install_binaries = rf"""
# Install CP2K binaries
COPY --from=build /opt/cp2k/exe/{arch}/cp2k.{version} \
                  /opt/cp2k/exe/{arch}/dumpdcd.{version} \
                  /opt/cp2k/exe/{arch}/graph.{version} \
                  /opt/cp2k/exe/{arch}/xyz2dcd.{version} \
                  /opt/cp2k/exe/{arch}/
"""
        run_tests = ""
        run_tests = r"""
# Create shortcut for regression test
RUN printf "echo Sorry, this container was built without test files" \
>/usr/local/bin/run_tests && chmod 755 /usr/local/bin/run_tests
"""
        tagname += "_no_tests"
    else:
        install_binaries = rf"""
# Install CP2K binaries
COPY --from=build /opt/cp2k/exe/{arch}/ /opt/cp2k/exe/{arch}/

# Install CP2K regression tests
COPY --from=build /opt/cp2k/tests/ /opt/cp2k/tests/
COPY --from=build /opt/cp2k/tools/regtesting/ /opt/cp2k/tools/regtesting/
COPY --from=build /opt/cp2k/src/grid/sample_tasks/ /opt/cp2k/src/grid/sample_tasks/
"""
        run_tests = rf"""
# Create shortcut for regression test
RUN printf "{do_regtest} {testopts} \$* {arch} {version}" \
>/usr/local/bin/run_tests && chmod 755 /usr/local/bin/run_tests
"""
        required_packages += " python3"

    return rf"""
# Usage: docker build -f ./Dockerfile.{name} -t cp2k/cp2k:{tagname} .

# Stage 1: build step
FROM {base_system} AS build

{cuda_environment}
# Install packages required for the CP2K toolchain build
RUN apt-get update -qq{permissive} && apt-get install -qq --no-install-recommends \
    {required_packages}{cuda_packages} \
    bzip2 ca-certificates git make patch pkg-config unzip wget zlib1g-dev

# Download CP2K
RUN git clone --recursive{branch} https://github.com/cp2k/cp2k.git /opt/cp2k

# Build CP2K toolchain for target CPU {target_cpu}
WORKDIR /opt/cp2k/tools/toolchain
RUN /bin/bash -c -o pipefail \
    "./install_cp2k_toolchain.sh -j {ncores} \
     --install-all \
     {with_gpu_line} \
     --target-cpu={target_cpu} \
     --with-cusolvermp=no \
     {with_gcc_line} \
     {with_mpi_line}"

# Build CP2K for target CPU {target_cpu}
WORKDIR /opt/cp2k
RUN /bin/bash -c -o pipefail \
    "cp ./tools/toolchain/install/arch/{arch}.{version} ./arch/; \
     source ./tools/toolchain/install/setup; \
     make -j {ncores} ARCH={arch} VERSION={version}{make_target}"

# Collect components for installation and remove symbolic links
RUN /bin/bash -c -o pipefail \
    "mkdir -p /toolchain/install /toolchain/scripts; \
     for libdir in \$(ldd ./exe/{arch}/cp2k.{version} | \
                      grep /opt/cp2k/tools/toolchain/install | \
                      awk '{{print \$3}}' | cut -d/ -f7 | \
                      sort | uniq) setup; do \
        cp -ar /opt/cp2k/tools/toolchain/install/\${{libdir}} /toolchain/install; \
     done; \
     cp /opt/cp2k/tools/toolchain/scripts/tool_kit.sh /toolchain/scripts; \
     unlink ./exe/{arch}/cp2k.{version.replace("smp", "opt")}; \
     unlink ./exe/{arch}/cp2k_shell.{version}"

# Stage 2: install step
FROM {base_system} AS install

# Install required packages
RUN apt-get update -qq{permissive} && apt-get install -qq --no-install-recommends \
    {required_packages} && rm -rf /var/lib/apt/lists/*
{install_binaries}
# Install CP2K database files
COPY --from=build /opt/cp2k/data/ /opt/cp2k/data/

# Install shared libraries required by the CP2K binaries
COPY --from=build /toolchain/ /opt/cp2k/tools/toolchain/

# Create links to CP2K binaries
RUN /bin/bash -c -o pipefail \
    "for binary in cp2k dumpdcd graph xyz2dcd; do \
        ln -sf /opt/cp2k/exe/{arch}/\${{binary}}.{version} \
               /usr/local/bin/\${{binary}}; \
     done; \
     ln -sf /opt/cp2k/exe/{arch}/cp2k.{version} \
            /usr/local/bin/cp2k_shell; \
     ln -sf /opt/cp2k/exe/{arch}/cp2k.{version} \
            /usr/local/bin/cp2k.{version.replace("smp", "opt")}"

# Create entrypoint script file
RUN printf "#!/bin/bash\n\
ulimit -c 0 -s unlimited\n\
{additional_exports}
export OMP_STACKSIZE={omp_stacksize}\n\
export PATH=/opt/cp2k/exe/{arch}:\${{PATH}}\n\
source /opt/cp2k/tools/toolchain/install/setup\n\
\"\$@\"" \
>/usr/local/bin/entrypoint.sh && chmod 755 /usr/local/bin/entrypoint.sh
{run_tests}
# Define entrypoint
WORKDIR /mnt
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["cp2k", "--help"]

# Label docker image
LABEL author="CP2K Developers" \
      cp2k_version="{release}" \
      dockerfile_generator_version="0.2"

# EOF
"""


# ------------------------------------------------------------------------------


class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()
        self.content.write("#\n")
        self.content.write("# This file was created by generate_docker_files.py\n")
        self.content.write("#")

    def __enter__(self) -> io.StringIO:
        return self.content

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        output_path = Path(__file__).parent / self.filename
        if self.check:
            assert output_path.read_text(encoding="utf8") == self.content.getvalue()
            print(f"File {output_path} is consistent with generator script")
        else:
            output_path.write_text(self.content.getvalue(), encoding="utf8")
            print(f"Wrote {output_path}")


# ------------------------------------------------------------------------------

main()

# EOF
