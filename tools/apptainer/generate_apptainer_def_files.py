#!/usr/bin/env python3

# Author: Matthias Krack (March 17, 2023)

from pathlib import Path
from typing import Any
import argparse
import io

# ------------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    distro = "ubuntu"
    distro_version = "22.04"
    arch = "local"
    version = "psmp"
    mpi_implementations = ["mpich"]
    omp_stacksize = "16M"
    testopts = "--workbasedir /mnt"
    version = "psmp"
    target_cpus = ["generic", "haswell", "broadwell", "skylake", "skylake-avx512"]
    for mpi_implementation in mpi_implementations:
        for target_cpu in target_cpus:
            with OutputFile(
                f"{mpi_implementation}_{target_cpu}_{version}.def", args.check
            ) as f:
                f.write(
                    write_definition_file(
                        distro=distro,
                        distro_version=distro_version,
                        arch=arch,
                        version=version,
                        mpi_implementation=mpi_implementation,
                        omp_stacksize=omp_stacksize,
                        target_cpu=target_cpu,
                        testopts=testopts,
                    )
                )


# ------------------------------------------------------------------------------


def write_definition_file(
    distro: str,
    distro_version: str,
    arch: str,
    version: str,
    mpi_implementation: str,
    omp_stacksize: str,
    target_cpu: str,
    testopts: str,
) -> str:
    return rf"""
# Usage: apptainer build -B $PWD:/mnt cp2k_{mpi_implementation}_{target_cpu}_{version}.sif {mpi_implementation}_{target_cpu}_{version}.def

Bootstrap: docker
From: {distro}:{distro_version}

%environment
 export CP2K_DATA_DIR=/opt/cp2k/data
 export OMP_STACKSIZE={omp_stacksize}

%post
 apt-get update -qq && apt-get install -qq --no-install-recommends \
  bzip2 ca-certificates g++ gcc gfortran git make pkg-config python3 unzip wget zlib1g-dev
 git clone --recursive https://github.com/cp2k/cp2k.git /opt/cp2k
 cd /opt/cp2k/tools/toolchain && ./install_cp2k_toolchain.sh \
  --install-all \
  --mpi-mode={mpi_implementation} \
  --target-cpu={target_cpu} \
  --with-gcc=system
 cd /opt/cp2k && cp ./tools/toolchain/install/arch/{arch}.{version} arch/
 bash -c -o pipefail \
 "source ./tools/toolchain/install/setup && \
  make -j ARCH={arch} VERSION={version} && \
  for binary in cp2k cp2k_shell dumpdcd graph xyz2dcd; do \
     ln -sf /opt/cp2k/exe/{arch}/\${{binary}}.{version} /usr/{arch}/bin/\${{binary}}; \
  done && \
  make -j ARCH={arch} VERSION={version} clean"
 cat ./tools/toolchain/install/setup >>$APPTAINER_ENVIRONMENT
 rm -rf ./tools/toolchain/build ./lib/{arch}/{version}/*.a ./.git

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
 Version v0.1
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
