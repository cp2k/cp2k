#!/usr/bin/env python3

# author: Ole Schuett

from pathlib import Path
from typing import Any, Optional
import argparse
import io

# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    for version in "sdbg", "ssmp", "pdbg", "psmp":
        with OutputFile(f"Dockerfile.test_{version}", args.check) as f:
            f.write(toolchain_full() + regtest(version))

        with OutputFile(f"Dockerfile.prod_{version}", args.check) as f:
            f.write(toolchain_full() + production(version))

    with OutputFile(f"Dockerfile.test_openmpi-psmp", args.check) as f:
        f.write(toolchain_full(mpi_mode="openmpi") + regtest("psmp"))

    with OutputFile(f"Dockerfile.test_fedora-psmp", args.check) as f:
        f.write(toolchain_full(base_image="fedora:33") + regtest("psmp"))

    with OutputFile(f"Dockerfile.test_minimal", args.check) as f:
        f.write(toolchain_full() + regtest("sdbg", "minimal"))

    for version in "sdbg", "pdbg":
        with OutputFile(f"Dockerfile.test_coverage-{version}", args.check) as f:
            f.write(toolchain_full() + coverage(version))

    for gcc_version in 7, 8, 9, 10:
        with OutputFile(f"Dockerfile.test_gcc{gcc_version}", args.check) as f:
            f.write(toolchain_ubuntu_nompi(gcc_version=gcc_version) + regtest("ssmp"))

    with OutputFile("Dockerfile.test_i386", args.check) as f:
        f.write(toolchain_ubuntu_nompi(base_image="i386/debian:11") + regtest("ssmp"))

    with OutputFile(f"Dockerfile.test_performance", args.check) as f:
        f.write(toolchain_full() + performance())

    for gpu_ver in "P100", "V100", "A100":
        with OutputFile(f"Dockerfile.test_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_cuda(gpu_ver=gpu_ver) + regtest("psmp", "local_cuda"))

        with OutputFile(f"Dockerfile.prod_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_cuda(gpu_ver=gpu_ver) + production("psmp", "local_cuda"))

        with OutputFile(f"Dockerfile.test_hip_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_hip_cuda(gpu_ver=gpu_ver) + regtest("psmp", "local_hip"))

        with OutputFile(f"Dockerfile.test_performance_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_cuda(gpu_ver=gpu_ver) + performance("local_cuda"))

    for gpu_ver in "Mi50", "Mi100":
        with OutputFile(f"Dockerfile.test_hip_rocm_{gpu_ver}", args.check) as f:
            f.write(toolchain_hip_rocm(gpu_ver=gpu_ver) + regtest("psmp", "local_hip"))

    with OutputFile(f"Dockerfile.test_conventions", args.check) as f:
        f.write(toolchain_full() + conventions())

    with OutputFile(f"Dockerfile.test_manual", args.check) as f:
        f.write(toolchain_full() + manual())

    with OutputFile(f"Dockerfile.test_precommit", args.check) as f:
        f.write(precommit())

    for name in "aiida", "ase", "gromacs", "i-pi":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(toolchain_full() + test_3rd_party(name))

    for name in "python", "doxygen":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(test_without_build(name))


# ======================================================================================
def regtest(version: str, arch: str = "local") -> str:
    return (
        install_cp2k(version=version, arch=arch)
        + fr"""
# Run regression tests.
ARG TESTOPTS
COPY ./tools/docker/scripts/test_regtest.sh ./
RUN /bin/bash -c " \
    TESTOPTS="${{TESTOPTS}}" \
    ./test_regtest.sh '{arch}' '{version}' |& tee report.log && \
    rm -rf regtesting"
"""
        + print_cached_report()
    )


# ======================================================================================
def performance(arch: str = "local") -> str:
    return (
        install_cp2k(version="psmp", arch=arch)
        + fr"""
# Run performance test for {arch}.
COPY ./benchmarks ./benchmarks
COPY ./tools/docker/scripts/test_performance.sh  \
     ./tools/docker/scripts/plot_performance.py  \
     ./
RUN ./test_performance.sh "{arch}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def coverage(version: str) -> str:
    return (
        install_cp2k(version=version, arch="local_coverage")
        + fr"""
# Run coverage test for {version}.
COPY ./tools/docker/scripts/test_coverage.sh .
RUN ./test_coverage.sh "{version}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def conventions() -> str:
    return (
        install_cp2k(version="psmp", arch="local_warn")
        + fr"""
# Run test for conventions.
COPY ./arch/Linux-x86-64-gfortran.dumpast ./arch/
COPY ./tools/conventions ./tools/conventions
WORKDIR ./tools/conventions
RUN /bin/bash -ec " \
    source /opt/cp2k-toolchain/install/setup && \
   ./test_conventions.sh |& tee report.log"
"""
        + print_cached_report()
    )


# ======================================================================================
def manual() -> str:
    return (
        install_cp2k(version="psmp", arch="local", revision=True)
        + fr"""
# Generate manual.
COPY ./tools/manual ./tools/manual
COPY ./tools/input_editing ./tools/input_editing
COPY ./tools/docker/scripts/test_manual.sh .
RUN ./test_manual.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def precommit() -> str:
    return (
        fr"""
FROM ubuntu:20.04

# Install dependencies.
WORKDIR /opt/cp2k-precommit
COPY ./tools/precommit/ /opt/cp2k-precommit/
RUN ./install_requirements.sh

# Install sources.
WORKDIR /opt/cp2k
COPY ./ ./

# Run precommit test.
RUN ./tools/docker/scripts/test_precommit.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def test_3rd_party(name: str) -> str:
    return (
        install_cp2k(version="sdbg", arch="local")
        + fr"""
# Run test for {name}.
COPY ./tools/docker/scripts/test_{name}.sh .
RUN ./test_{name}.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def test_without_build(name: str) -> str:
    return (
        fr"""
FROM ubuntu:20.04

# Install dependencies.
WORKDIR /opt/cp2k
COPY ./tools/docker/scripts/install_{name}.sh .
RUN ./install_{name}.sh

# Install sources.
ARG GIT_COMMIT_SHA
COPY ./src ./src
COPY ./exts ./exts
COPY ./data ./data
COPY ./tools ./tools
COPY ./Makefile .
RUN bash -c "if [ -n "${{GIT_COMMIT_SHA}}" ] ; then echo "git:\${{GIT_COMMIT_SHA::7}}" > REVISION; fi"

# Run test for {name}.
COPY ./tools/docker/scripts/test_{name}.sh .
RUN ./test_{name}.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def print_cached_report() -> str:
    return r"""
# Output the report if the image is old and was therefore pulled from the build cache.
CMD cat $(find ./report.log -mmin +10) | sed '/^Summary:/ s/$/ (cached)/'

#EOF
"""


# ======================================================================================
def production(version: str, arch: str = "local") -> str:
    return (
        install_cp2k(version=version, arch=arch, revision=True, prod=True)
        + fr"""
# Run regression tests.
ARG TESTOPTS
RUN /bin/bash -c " \
    source /opt/cp2k-toolchain/install/setup && \
    ./tools/regtesting/do_regtest.py '{arch}' '{version}' "${{TESTOPTS}}" |& tee regtests.log && \
    rm -rf regtesting"

# Setup entry point for production.
COPY ./tools/docker/scripts/prod_entrypoint.sh ./
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/prod_entrypoint.sh", "{arch}", "{version}"]
CMD ["cp2k", "--help"]

#EOF
"""
    )


# ======================================================================================
def install_cp2k(
    version: str, arch: str, revision: bool = False, prod: bool = False
) -> str:
    input_lines = []
    run_lines = []

    if revision:
        input_lines.append("ARG GIT_COMMIT_SHA")
        run_lines.append(
            'if [ -n "${GIT_COMMIT_SHA}" ] ; then'
            ' echo "git:\${GIT_COMMIT_SHA::7}" > REVISION; fi'
        )

    input_lines.append("COPY ./Makefile .")
    input_lines.append("COPY ./src ./src")
    input_lines.append("COPY ./exts ./exts")
    input_lines.append("COPY ./tools/build_utils ./tools/build_utils")

    if arch.startswith("local"):
        arch_file = f"/opt/cp2k-toolchain/install/arch/{arch}.{version}"
        run_lines.append("mkdir -p arch")
        run_lines.append(f"ln -vs {arch_file} ./arch/")
    else:
        input_lines.append(f"COPY ./arch/{arch}.{version} /opt/cp2k/arch/")

    run_lines.append("echo 'Compiling cp2k...'")
    run_lines.append("source /opt/cp2k-toolchain/install/setup")

    if prod:
        run_lines.append(f"make -j ARCH={arch} VERSION={version}")
        run_lines.append(f"rm -rf lib obj exe/{arch}/libcp2k_unittest.{version}")
    else:
        run_lines.append(
            f"( make -j ARCH={arch} VERSION={version} &> /dev/null || true )"
        )

    input_block = "\n".join(input_lines)
    run_block = " && \\\n    ".join(run_lines)

    return fr"""
# Install CP2K using {arch}.{version}.
WORKDIR /opt/cp2k
{input_block}
RUN /bin/bash -c " \
    {run_block}"
COPY ./data ./data
COPY ./tests ./tests
COPY ./tools/regtesting ./tools/regtesting
"""


# ======================================================================================
def toolchain_full(base_image: str = "ubuntu:20.04", mpi_mode: str = "mpich") -> str:
    return f"\nFROM {base_image}\n\n" + install_toolchain(
        base_image=base_image, install_all=None, mpi_mode=mpi_mode
    )


# ======================================================================================
def toolchain_ubuntu_nompi(
    base_image: str = "ubuntu:20.04", gcc_version: int = 10
) -> str:
    return fr"""
FROM {base_image}

# Install Ubuntu packages.
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

""" + install_toolchain(
        base_image="ubuntu",
        mpi_mode="no",
        with_gcc="system",
        with_cmake="system",
        with_fftw="system",
        with_openblas="system",
        with_gsl="system",
        with_hdf5="system",
        with_libxc="install",
        with_libxsmm="install",
        with_libint="install",
    )


# ======================================================================================
def toolchain_cuda(gpu_ver: str) -> str:
    return fr"""
FROM nvidia/cuda:11.3.1-devel-ubuntu20.04

# Setup CUDA environment.
ENV CUDA_PATH /usr/local/cuda
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64

# Install Ubuntu packages.
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    gfortran                                                          \
    mpich                                                             \
    libmpich-dev                                                      \
   && rm -rf /var/lib/apt/lists/*

""" + install_toolchain(
        base_image="ubuntu", mpi_mode="mpich", enable_cuda="yes", gpu_ver=gpu_ver
    )


# ======================================================================================
def toolchain_hip_cuda(gpu_ver: str) -> str:
    return fr"""
FROM nvidia/cuda:11.3.1-devel-ubuntu20.04

# Setup CUDA environment.
ENV CUDA_PATH /usr/local/cuda
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64
ENV HIP_PLATFORM nvidia
ENV ROCM_VER 4.5.2
ENV HIP_DIR /opt/HIP-rocm-4.5.2
ENV HIPAMD_DIR /opt/hipamd-rocm-4.5.2

RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true \
    && apt-get update -qq && apt-get install -qq --no-install-recommends \
    ca-certificates \
    build-essential \
    cmake \
    git \
    gfortran \
    mpich \
    libmpich-dev \
    wget \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Install HIP from source because the hip-nvcc package drags in 10GB of unnecessary dependencies.
WORKDIR /opt

RUN wget -q https://github.com/Kitware/CMake/releases/download/v3.20.6/cmake-3.20.6-Linux-x86_64.sh \
    && echo "4772100c2578927eed5aa9e1a80694c0d64410448c0fda73d31b0eae18645784  cmake-3.20.6-Linux-x86_64.sh" | sha256sum --check \
    && sh cmake-3.20.6-Linux-x86_64.sh --prefix=/usr/local --skip-license \
    && rm -f cmake-3.20.6-Linux-x86_64.sh \
    && cmake --version

RUN wget -q https://github.com/ROCm-Developer-Tools/HIP/archive/refs/tags/rocm-${{ROCM_VER}}.tar.gz -O HIP-rocm-${{ROCM_VER}}.tar.gz\
    && echo "c2113dc3c421b8084cd507d91b6fbc0170765a464b71fb0d96bb875df368f160  HIP-rocm-${{ROCM_VER}}.tar.gz" |  sha256sum --check \
    && tar -xzf HIP-rocm-*.tar.gz \
    && wget -q https://github.com/ROCm-Developer-Tools/hipamd/archive/refs/tags/rocm-${{ROCM_VER}}.tar.gz -O hipamd-rocm-${{ROCM_VER}}.tar.gz \
    && echo "b6f35b1a1d0c466b5af28e26baf646ae63267eccc4852204db1e0c7222a39ce2  hipamd-rocm-${{ROCM_VER}}.tar.gz" | sha256sum --check \
    && tar -xzf hipamd-rocm-*.tar.gz \
    && wget -q https://github.com/ROCmSoftwarePlatform/hipBLAS/archive/refs/tags/rocm-${{ROCM_VER}}.tar.gz -O hipBLAS-rocm-${{ROCM_VER}}.tar.gz \
    && echo "82dd82a41bbadbb2a91a2a44a5d8e0d2e4f36d3078286ed4db3549b1fb6d6978  hipBLAS-rocm-${{ROCM_VER}}.tar.gz" | sha256sum --check \
    && tar -xzf hipBLAS-rocm-*.tar.gz \
    && wget -q https://github.com/ROCmSoftwarePlatform/hipFFT/archive/refs/tags/rocm-${{ROCM_VER}}.tar.gz -O hipFFT-rocm-${{ROCM_VER}}.tar.gz \
    && echo "32ba6a5f50cfede3777a43794371ffb1363302131d8a0382d96df90ed7bc911a  hipFFT-rocm-${{ROCM_VER}}.tar.gz" | sha256sum --check \
    && tar -xzf hipFFT-rocm-*.tar.gz

RUN cd ${{HIPAMD_DIR}} \
    && mkdir -p build \
    && cd  build \
    && mkdir /opt/rocm-${{ROCM_VER}} \
    && cmake -DHIP_COMMON_DIR=${{HIP_DIR}} -DHIP_PLATFORM=nvidia -DCMAKE_INSTALL_PREFIX=/opt/rocm-${{ROCM_VER}}/hip .. > /dev/null 2>&1 \
    && make -j > /dev/null 2>&1 \
    && make install > /dev/null 2>&1 \
    && cd ../..

# Install hipBLAS from source.
RUN cd hipBLAS-rocm-* \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/rocm-${{ROCM_VER}} -DUSE_CUDA=YES -DCMAKE_MODULE_PATH=/opt/rocm-${{ROCM_VER}} -DCMAKE_MODULE_PATH=/opt/rocm-${{ROCM_VER}}/hip/cmake .. > /dev/null 2>&1 \
    && make -j > /dev/null 2>&1 \
    && make install > /dev/null 2>&1 \
    && cd .. \
    && rm -rf hipBLAS-rocm-*

ENV CPATH ${{CPATH}}:/opt/rocm-${{ROCM_VER}}/hip/include
# Install hipFFT from source.
RUN cd hipFFT-rocm-* \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/opt/rocm-${{ROCM_VER}} -DBUILD_WITH_LIB=CUDA .. > /dev/null 2>&1 \
    && make -j > /dev/null 2>&1 \
    && make install > /dev/null 2>&1 \
    && rm -rf hipFFT*

# This is the alternative installation path via Ubuntu packages.
## https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation-Guide.html#ubuntu
## https://rocmdocs.amd.com/en/latest/Installation_Guide/HIP-Installation.html#nvidia-platform
#RUN apt-key adv --fetch-keys https://repo.radeon.com/rocm/rocm.gpg.key
#RUN echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/debian/ ubuntu main' > /etc/apt/sources.list.d/rocm.list
#RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true \
#    && apt-get update -qq \
#    && apt-get install --yes --no-install-recommends hip-nvcc hipblas \
#    && rm -rf /var/lib/apt/lists/*

# Setup HIP environment.
ENV ROCM_PATH /opt/rocm-${{ROCM_VER}}
ENV PATH ${{PATH}}:${{ROCM_PATH}}/hip/bin
ENV LD_LIBRARY_PATH ${{LD_LIBRARY_PATH}}:${{ROCM_PATH}}/hip/lib:${{ROCM_PATH}}/hipblas/lib:${{ROCM_PATH}}/hipfft/lib:${{ROCM_PATH}}/hipfft/lib64:${{ROCM_PATH}}/hipblas/lib64
ENV HIP_PLATFORM nvidia
RUN hipconfig

""" + install_toolchain(
        base_image="ubuntu", mpi_mode="mpich", enable_hip="yes", gpu_ver=gpu_ver
    )


# ======================================================================================
def toolchain_hip_rocm(gpu_ver: str) -> str:
    return fr"""
FROM rocm/dev-ubuntu-20.04:4.5.2-complete

# Install some Ubuntu packages.
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    hipblas                                                           \
    gfortran                                                          \
    mpich                                                             \
    libmpich-dev                                                      \
   && rm -rf /var/lib/apt/lists/*

# Setup HIP environment.
ENV ROCM_PATH /opt/rocm
ENV PATH ${{PATH}}:${{ROCM_PATH}}/bin
ENV LD_LIBRARY_PATH ${{LD_LIBRARY_PATH}}:${{ROCM_PATH}}/lib
ENV HIP_PLATFORM amd
RUN hipconfig

""" + install_toolchain(
        base_image="ubuntu", mpi_mode="mpich", enable_hip="yes", gpu_ver=gpu_ver
    )


# ======================================================================================
def install_toolchain(base_image: str, **kwargs: Optional[str]) -> str:
    install_args = []
    for k, v in kwargs.items():
        k = k.replace("_", "-")
        if v is not None:
            install_args.append(f"    --{k}={v} \\")
        else:
            install_args.append(f"    --{k} \\")
    install_args_str = "\n".join(install_args)

    return fr"""
# Install requirements for the toolchain.
WORKDIR /opt/cp2k-toolchain
COPY ./tools/toolchain/install_requirements*.sh ./
RUN ./install_requirements.sh {base_image}

# Install the toolchain.
RUN mkdir scripts
COPY ./tools/toolchain/scripts/VERSION \
     ./tools/toolchain/scripts/parse_if.py \
     ./tools/toolchain/scripts/tool_kit.sh \
     ./tools/toolchain/scripts/common_vars.sh \
     ./tools/toolchain/scripts/signal_trap.sh \
     ./tools/toolchain/scripts/get_openblas_arch.sh \
     ./scripts/
COPY ./tools/toolchain/install_cp2k_toolchain.sh .
RUN ./install_cp2k_toolchain.sh \
{install_args_str}
    --dry-run

# Dry-run leaves behind config files for the followup install scripts.
# This breaks up the lengthy installation into smaller docker build steps.
COPY ./tools/toolchain/scripts/stage0/ ./scripts/stage0/
RUN  ./scripts/stage0/install_stage0.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage1/ ./scripts/stage1/
RUN  ./scripts/stage1/install_stage1.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage2/ ./scripts/stage2/
RUN  ./scripts/stage2/install_stage2.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage3/ ./scripts/stage3/
RUN  ./scripts/stage3/install_stage3.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage4/ ./scripts/stage4/
RUN  ./scripts/stage4/install_stage4.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage5/ ./scripts/stage5/
RUN  ./scripts/stage5/install_stage5.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage6/ ./scripts/stage6/
RUN  ./scripts/stage6/install_stage6.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage7/ ./scripts/stage7/
RUN  ./scripts/stage7/install_stage7.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/stage8/ ./scripts/stage8/
RUN  ./scripts/stage8/install_stage8.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/arch_base.tmpl \
     ./tools/toolchain/scripts/generate_arch_files.sh \
     ./scripts/
RUN ./scripts/generate_arch_files.sh && rm -rf ./build
""".lstrip()


# ======================================================================================
class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()
        self.content.write(f"#\n")
        self.content.write(f"# This file was created by generate_dockerfiles.py.\n")
        self.content.write(f"# Usage: docker build -f ./{filename} ../../\n")
        self.content.write(f"#\n")

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
