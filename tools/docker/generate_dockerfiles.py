#!/usr/bin/env python3

# author: Ole Schuett

from pathlib import Path
from typing import Any
import argparse
import io

# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    for gcc_version in 7, 8, 9, 10:
        with OutputFile(f"Dockerfile.test_gcc{gcc_version}", args.check) as f:
            f.write(toolchain_ubuntu_nompi(gcc_version=gcc_version) + regtest_ssmp())

    with OutputFile("Dockerfile.test_i386", args.check) as f:
        f.write(toolchain_ubuntu_nompi(base_image="i386/debian:11") + regtest_ssmp())


# ======================================================================================
def regtest_ssmp() -> str:
    return fr"""
WORKDIR /workspace

COPY ./tools/docker/scripts/install_basics.sh .
RUN ./install_basics.sh

COPY ./tools/docker/scripts/install_regtest.sh .
RUN ./install_regtest.sh local ssmp

COPY ./tools/docker/scripts/ci_entrypoint.sh ./tools/docker/scripts/test_regtest.sh ./
CMD ["./ci_entrypoint.sh", "./test_regtest.sh", "local", "ssmp"]

#EOF
""".lstrip()


# ======================================================================================
def toolchain_ubuntu_nompi(
    base_image: str = "ubuntu:20.04", gcc_version: int = 10, libint_lmax: int = 4
) -> str:
    return fr"""
FROM {base_image}
USER root

# author: Ole Schuett

# Installs lean toolchain without MPI and relying mostly on Ubuntu packages.

# Install Ubuntu packages.
COPY ./tools/toolchain/install_requirements_ubuntu.sh .
RUN ./install_requirements_ubuntu.sh

# Install some more Ubuntu packages.
RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
    apt-get update -qq && apt-get install -qq --no-install-recommends \
    cmake \
    gcc-{gcc_version} \
    g++-{gcc_version} \
    gfortran-{gcc_version} \
    fftw3-dev \
    libopenblas-dev \
    libgsl-dev \
    libhdf5-dev \
   && rm -rf /var/lib/apt/lists/*

# Create links.
RUN ln -sf gcc-{gcc_version}      /usr/bin/gcc  && \
    ln -sf g++-{gcc_version}      /usr/bin/g++  && \
    ln -sf gfortran-{gcc_version} /usr/bin/gfortran

# Build toolchain.
WORKDIR /opt/cp2k-toolchain
COPY ./tools/toolchain/scripts ./scripts
COPY ./tools/toolchain/install_cp2k_toolchain.sh .
RUN ./install_cp2k_toolchain.sh \
    --mpi-mode=no \
    --with-gcc=system \
    --with-cmake=system \
    --with-fftw=system \
    --with-openblas=system \
    --with-gsl=system \
    --with-hdf5=system \
    --with-libxc=install \
    --with-libxsmm=install \
    --with-libint=install \
     --with-libxc=install \
    --with-libxsmm=install \
    --with-libint=install \
    --libint-lmax={libint_lmax} \
    && rm -rf ./build
""".lstrip()


# ======================================================================================
class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()

    def __enter__(self) -> io.StringIO:
        return self.content

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        output_path = Path(__file__).parent / self.filename
        if self.check:
            assert output_path.read_text(encoding="utf8") == self.content.getvalue()
            print(f"File {output_path} is consisted with generator script.")
        else:
            output_path.write_text(self.content.getvalue(), encoding="utf8")
            print(f"Wrote {output_path}")


# ======================================================================================
main()

# EOF
