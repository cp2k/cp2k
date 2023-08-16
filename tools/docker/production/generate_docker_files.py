#!/usr/bin/env python3

# Author: Matthias Krack (August 15, 2023)

from pathlib import Path
from typing import Any
import argparse
import io
import os

# ------------------------------------------------------------------------------

cp2k_release_list = ["master", "2023.2"]  # append new releases to list
mpi_implementation_list = ["mpich", "openmpi"]
mpich_device_list = ["ch3", "ch4", "ch4:ucx"]
target_cpu_list = ["generic", "haswell", "skylake-avx512", "native", "znver2", "znver3"]

# ------------------------------------------------------------------------------


def main() -> None:
    mpi_choices = ["all"] + mpi_implementation_list
    release_choices = ["all"] + cp2k_release_list
    target_cpu_choices = ["all"] + target_cpu_list

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check",
        action="store_true",
        dest="check",
        help="Check consistency with generator script",
    )
    parser.add_argument(
        "--install-tests",
        action="store_true",
        dest="install_tests",
        help="Install tests and enable regression testing for the final container",
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
        "--mpich-device",
        choices=mpich_device_list,
        default=mpich_device_list[-1],
        dest="mpich_device",
        help=f"Select the MPICH device (default is {mpich_device_list[-1]})",
        type=str,
    )
    parser.add_argument(
        "-j",
        "--ncores",
        default=8,
        dest="ncores",
        help=(
            "Select the number CPU cores used for building the container "
            "(default is 8)"
        ),
        type=check_ncores,
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
            "Specify the CP2K release for which docker files are generated "
            f"(default is {target_cpu_choices[0]})"
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
    args = parser.parse_args()

    distro = "ubuntu"
    distro_version = "22.04"
    arch = "local"
    version = "psmp"
    mpich_device = args.mpich_device
    ncores = args.ncores
    omp_stacksize = "16M"
    install_tests = args.install_tests
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

            for target_cpu in target_cpu_list:
                if args.target_cpu != "all" and args.target_cpu != target_cpu:
                    continue
                name = f"{release}_{mpi_implementation}_{target_cpu}_{version}"
                with OutputFile(f"Dockerfile.{name}", args.check) as f:
                    f.write(
                        write_definition_file(
                            name=name,
                            release=release,
                            distro=distro,
                            distro_version=distro_version,
                            arch=arch,
                            version=version,
                            mpi_implementation=mpi_implementation,
                            mpich_device=mpich_device,
                            ncores=ncores,
                            omp_stacksize=omp_stacksize,
                            install_tests=install_tests,
                            target_cpu=target_cpu,
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
    name: str,
    release: str,
    distro: str,
    distro_version: str,
    arch: str,
    version: str,
    mpi_implementation: str,
    mpich_device: str,
    ncores: int,
    omp_stacksize: str,
    install_tests: bool,
    target_cpu: str,
    test_build: bool,
) -> str:
    if release == "master":
        branch = ""
        tagname = name.replace("master", r"master$(date +%Y%m%d)")
    else:
        branch = f" -b support/v{release}"
        tagname = name

    if test_build:
        make_target = " test"
    else:
        make_target = ""

    additional_exports = "\\"
    with_mpi_line = ""

    # Required packages for final container
    required_packages = "g++ gcc gfortran"

    # Default options for the regression tests
    testopts = f"--maxtasks {ncores}" " --workbasedir /mnt"

    if mpi_implementation == "openmpi":
        additional_exports = """\
export OMPI_ALLOW_RUN_AS_ROOT=1\\n\\
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\\n\\
export OMPI_MCA_btl_vader_single_copy_mechanism=none\\n\\\
"""
        with_mpi_line = f"--with-{mpi_implementation}=install"
        required_packages += " openssh-client"
        testopts += '--mpiexec "mpiexec --bind-to none" ' + testopts
    elif mpi_implementation == "mpich":
        with_mpi_line = (
            f"--with-{mpi_implementation}=install "
            f"--with-{mpi_implementation}-device={mpich_device}"
        )

    if install_tests:
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
RUN printf "/opt/cp2k/tools/regtesting/do_regtest.py {testopts} \$* {arch} {version}" \
>/usr/local/bin/run_tests && chmod 755 /usr/local/bin/run_tests
"""
        required_packages += " python3"
    else:
        install_binaries = rf"""
# Install CP2K binaries
COPY --from=build /opt/cp2k/exe/{arch}/cp2k.{version} \
                  /opt/cp2k/exe/{arch}/dumpdcd.{version} \
                  /opt/cp2k/exe/{arch}/graph.{version} \
                  /opt/cp2k/exe/{arch}/xyz2dcd.{version} \
                  /opt/cp2k/exe/{arch}/
"""
        run_tests = ""

    return rf"""
# Usage: docker build -f ./Dockerfile.{name} -t cp2k/cp2k:{tagname} .

# Stage 1: build step
FROM {distro}:{distro_version} AS build

# Install packages required for the CP2K toolchain build
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    bzip2 ca-certificates g++ gcc gfortran git make openssh-client patch \
    pkg-config python3 unzip wget zlib1g-dev

# Download CP2K
RUN git clone --recursive{branch} https://github.com/cp2k/cp2k.git /opt/cp2k

# Build CP2K toolchain for target CPU {target_cpu}
WORKDIR /opt/cp2k/tools/toolchain
RUN /bin/bash -c -o pipefail \
    "./install_cp2k_toolchain.sh \
     --install-all \
     --target-cpu={target_cpu} \
     --with-gcc=system \
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
FROM ubuntu:22.04 AS install

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
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
export OMP_STACKSIZE={omp_stacksize}\n\
export PATH=/opt/cp2k/exe/{arch}:\${{PATH}}\n\
{additional_exports}
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
      dockerfile_generator_version="0.1"

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
