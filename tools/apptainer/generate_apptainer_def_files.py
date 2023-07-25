#!/usr/bin/env python3

# Author: Matthias Krack (March 24, 2023)

from pathlib import Path
from typing import Any
import argparse
import io
import os

# ------------------------------------------------------------------------------

cp2k_release_list = ["master", "2023.1"]  # append new releases to list
mpi_implementation_list = ["mpich", "openmpi"]
mpich_device_list = ["ch3", "ch4", "ch4:ucx"]
target_cpu_list = ["generic", "haswell", "skylake-avx512", "native", "znver2", "znver3"]


def main() -> None:
    mpi_choices = ["all"] + mpi_implementation_list
    release_choices = ["all"] + cp2k_release_list
    target_cpu_choices = ["all"] + target_cpu_list

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check", action="store_true", help="Check consistency with generator script"
    )
    parser.add_argument(
        "--mpi",
        choices=mpi_choices,
        default=mpi_choices[0],
        dest="mpi_implementation",
        help=f"Select a MPI implementation (default is to generate apptainer definition files for {mpi_choices[0]})",
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
        help="Select the number CPU cores used for building the container (default is 8)",
        type=check_ncores,
    )
    parser.add_argument(
        "--release",
        choices=release_choices,
        default=release_choices[0],
        dest="release",
        help=f"Specify the CP2K release for which definition files are generated (default is {release_choices[0]})",
        type=str,
    )
    parser.add_argument(
        "--target-cpu",
        choices=target_cpu_choices,
        default=target_cpu_choices[0],
        dest="target_cpu",
        help=f"Specify the CP2K release for which definition files are generated (default is {target_cpu_choices[0]})",
        type=str,
    )
    args = parser.parse_args()

    package = "cp2k"
    distro = "ubuntu"
    distro_version = "22.04"
    arch = "local"
    version = "psmp"
    mpich_device = args.mpich_device
    ncores = args.ncores
    omp_stacksize = "16M"
    testopts = f"--keepalive --maxtasks {ncores} --skipdir UNIT/libcp2k_unittest --workbasedir /mnt"
    version = "psmp"

    if ncores > os.cpu_count():
        print(
            f"WARNING: More CPU cores requested for build than available ({ncores} > {os.cpu_count()})"
        )

    for release in cp2k_release_list:
        if args.release != "all" and args.release != release:
            continue

        if release == "master":
            prefix = f"{package}"
        else:
            prefix = f"{package}-{release}"

        for mpi_implementation in mpi_implementation_list:
            if (
                args.mpi_implementation != "all"
                and args.mpi_implementation != mpi_implementation
            ):
                continue

            for target_cpu in target_cpu_list:
                if args.target_cpu != "all" and args.target_cpu != target_cpu:
                    continue
                if target_cpu == "haswell" and release == "2023.1":
                    continue  # skip because of a toolchain bug in v2023.1
                name = f"{mpi_implementation}_{target_cpu}_{version}"
                defname = f"{release}_{name}"
                sifname = f"{prefix}_{name}"
                with OutputFile(f"{defname}.def", args.check) as f:
                    f.write(
                        write_definition_file(
                            defname=defname,
                            sifname=sifname,
                            release=release,
                            distro=distro,
                            distro_version=distro_version,
                            arch=arch,
                            version=version,
                            mpi_implementation=mpi_implementation,
                            mpich_device=mpich_device,
                            ncores=ncores,
                            omp_stacksize=omp_stacksize,
                            target_cpu=target_cpu,
                            testopts=testopts,
                        )
                    )


# ------------------------------------------------------------------------------


def check_ncores(value: str) -> int:
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError("%s is an invalid number of CPU cores" % value)
    return ivalue


# ------------------------------------------------------------------------------


def write_definition_file(
    defname: str,
    sifname: str,
    release: str,
    distro: str,
    distro_version: str,
    arch: str,
    version: str,
    mpi_implementation: str,
    mpich_device: str,
    ncores: int,
    omp_stacksize: str,
    target_cpu: str,
    testopts: str,
) -> str:
    if release == "master":
        branch = ""
    else:
        branch = f" -b support/v{release}"

    additional_exports = ""
    with_mpi_line = ""

    if mpi_implementation == "openmpi":
        additional_exports = (
            "export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\n"
        )
        with_mpi_line = f"--with-{mpi_implementation}=install"
        testopts = f'--mpiexec "mpiexec --bind-to none" ' + testopts
    elif mpi_implementation == "mpich":
        with_mpi_line = f"--with-{mpi_implementation}=install --with-{mpi_implementation}-device={mpich_device}"

    return rf"""
# Usage: apptainer build -B $PWD:/mnt {sifname}.sif {defname}.def | tee {sifname}.log

Bootstrap: docker
From: {distro}:{distro_version}

%environment
 export CP2K_DATA_DIR=/opt/cp2k/data
 export OMP_STACKSIZE={omp_stacksize}
 {additional_exports}
%post
 apt-get update -qq && apt-get install -qq --no-install-recommends \
  bzip2 ca-certificates g++ gcc gfortran git make openssh-client patch pkg-config python3 unzip wget zlib1g-dev
 git clone --recursive{branch} https://github.com/cp2k/cp2k.git /opt/cp2k
 cd /opt/cp2k/tools/toolchain && ./install_cp2k_toolchain.sh -j {ncores} \
  --install-all \
  --target-cpu={target_cpu} \
  --with-gcc=system \
  {with_mpi_line}
 cd /opt/cp2k && cp ./tools/toolchain/install/arch/{arch}.{version} arch/
 bash -c -o pipefail \
 "source ./tools/toolchain/install/setup && \
  make -j {ncores} ARCH={arch} VERSION={version} && \
  for binary in cp2k cp2k_shell dumpdcd graph xyz2dcd; do \
     ln -sf /opt/cp2k/exe/{arch}/\${{binary}}.{version} /usr/{arch}/bin/\${{binary}}; \
  done && \
  make -j {ncores} ARCH={arch} VERSION={version} clean"
 cat ./tools/toolchain/install/setup >>$APPTAINER_ENVIRONMENT
 rm -rf ./tools/toolchain/build ./lib/{arch}/{version}/*.a ./exe/local/libcp2k_unittest.psmp ./.git

%test
 #!/bin/bash
 export OMP_STACKSIZE={omp_stacksize}
 ulimit -c 0 -s unlimited
 /opt/cp2k/tools/regtesting/do_regtest.py {testopts} {arch} {version}

%runscript
 #!/bin/bash
 ulimit -c 0 -s unlimited
 "$@"

%labels
 Author CP2K developer team
 Version v0.2
"""


# ------------------------------------------------------------------------------


class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()
        self.content.write(f"#\n")
        self.content.write(
            f"# This file was created by generate_apptainer_def_file.py\n"
        )
        self.content.write(f"#")

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
