#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import io
from pathlib import Path
from typing import Any


# ======================================================================================
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    for version in "sdbg", "ssmp", "pdbg", "psmp":
        with OutputFile(f"Dockerfile.test_{version}", args.check) as f:
            mpi_mode = "mpich" if version.startswith("p") else "no"
            f.write(install_deps_toolchain(mpi_mode=mpi_mode, with_dbcsr=""))
            f.write(regtest_cmake("toolchain_all", version))

    with OutputFile(f"Dockerfile.test_generic_psmp", args.check) as f:
        f.write(install_deps_toolchain(target_cpu="generic", with_dbcsr=""))
        f.write(regtest_cmake("toolchain_generic", "psmp"))

    with OutputFile(f"Dockerfile.test_openmpi-psmp", args.check) as f:
        f.write(install_deps_toolchain(mpi_mode="openmpi", with_dbcsr=""))
        f.write(regtest_cmake("toolchain_all", "psmp"))

    with OutputFile(f"Dockerfile.test_fedora-psmp", args.check) as f:
        f.write(install_deps_toolchain(base_image="fedora:41", with_dbcsr=""))
        f.write(regtest_cmake("toolchain_all", "psmp"))

    for ver in "ssmp", "psmp":
        with OutputFile(f"Dockerfile.test_intel-{ver}", args.check) as f:
            base_image = "intel/hpckit:2024.2.1-0-devel-ubuntu22.04"
            f.write(install_deps_toolchain_intel(base_image=base_image, with_ifx="no"))
            f.write(regtest(ver, intel=True, testopts="--mpiexec mpiexec"))
        with OutputFile(f"Dockerfile.test_intel-oneapi-hpckit-{ver}", args.check) as f:
            base_image = "intel/oneapi-hpckit:2025.1.3-0-devel-ubuntu24.04"
            f.write(install_deps_toolchain_intel(base_image=base_image, with_ifx="yes"))
            f.write(regtest(ver, intel=True, testopts="--mpiexec mpiexec"))

    with OutputFile(f"Dockerfile.test_nvhpc", args.check) as f:
        f.write(install_deps_toolchain_nvhpc())

    with OutputFile(f"Dockerfile.test_minimal", args.check) as f:
        f.write(install_deps_ubuntu())
        f.write(install_dbcsr("minimal", "ssmp"))
        f.write(regtest_cmake("minimal", "ssmp"))

    with OutputFile(f"Dockerfile.test_spack", args.check) as f:
        f.write(install_deps_spack("psmp"))
        f.write(regtest_cmake("spack_all", "psmp"))

    with OutputFile(f"Dockerfile.test_asan-psmp", args.check) as f:
        f.write(install_deps_toolchain(with_dbcsr=""))
        f.write(regtest_cmake("toolchain_asan", "psmp"))

    for version in "sdbg", "pdbg":
        with OutputFile(f"Dockerfile.test_coverage-{version}", args.check) as f:
            f.write(install_deps_toolchain())
            f.write(coverage(version))

    for gcc_version in 8, 9, 10, 11, 12, 13, 14:
        with OutputFile(f"Dockerfile.test_gcc{gcc_version}", args.check) as f:
            if gcc_version > 8:
                f.write(install_deps_ubuntu(gcc_version=gcc_version))
                f.write(install_dbcsr("ubuntu", "ssmp"))
                f.write(regtest_cmake("ubuntu", "ssmp"))
            else:
                f.write(install_deps_ubuntu2004(gcc_version=gcc_version))
                # Have to use Makefile because Ubuntu:20.04 ships with CMake 3.16.3.
                # Skip some tests due to bug in LDA_C_PMGB06 functional in libxc <5.2.0.
                f.write(regtest("ssmp", testopts="--skipdir=QS/regtest-rs-dhft"))

    with OutputFile("Dockerfile.test_arm64-psmp", args.check) as f:
        base_img = "arm64v8/ubuntu:24.04"
        f.write(install_deps_toolchain(base_img, with_libtorch="no", with_deepmd="no"))
        f.write(regtest("psmp"))

    with OutputFile(f"Dockerfile.test_performance", args.check) as f:
        f.write(install_deps_toolchain())
        f.write(performance())

    for gpu_ver in "P100", "V100", "A100":
        with OutputFile(f"Dockerfile.test_cuda_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_cuda(gpu_ver=gpu_ver))
            f.write(regtest("psmp", "local_cuda"))

        with OutputFile(f"Dockerfile.test_hip_cuda_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_hip_cuda(gpu_ver=gpu_ver))
            f.write(regtest("psmp", "local_hip"))

        with OutputFile(f"Dockerfile.test_performance_cuda_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_cuda(gpu_ver=gpu_ver))
            f.write(performance("local_cuda"))

    for gpu_ver in "Mi50", "Mi100":
        with OutputFile(f"Dockerfile.test_hip_rocm_{gpu_ver}", args.check) as f:
            # ROCm containers require --device, which is not available for docker build.
            # https://rocmdocs.amd.com/en/latest/ROCm_Virtualization_Containers/ROCm-Virtualization-&-Containers.html#docker-hub
            f.write(install_deps_toolchain_hip_rocm(gpu_ver=gpu_ver))
            f.write(regtest_postponed("psmp", "local_hip"))

        with OutputFile(f"Dockerfile.build_hip_rocm_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_hip_rocm(gpu_ver=gpu_ver))
            f.write(build("psmp", "local_hip"))

    with OutputFile(f"Dockerfile.test_conventions", args.check) as f:
        f.write(install_deps_toolchain())
        f.write(conventions())

    with OutputFile(f"Dockerfile.test_manual", args.check) as f:
        f.write(install_deps_toolchain())
        f.write(manual())

    with OutputFile(f"Dockerfile.test_precommit", args.check) as f:
        f.write(precommit())

    for name in "ase", "aiida", "i-pi", "phonopy", "gromacs":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(install_deps_toolchain(mpi_mode="no", with_dbcsr=""))
            f.write(test_3rd_party(name))

    for name in "misc", "doxygen":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(test_without_build(name))


# ======================================================================================
def regtest(
    version: str, arch: str = "local", testopts: str = "", intel: bool = False
) -> str:
    return (
        install_cp2k(version=version, arch=arch, intel=intel)
        + rf"""
# Run regression tests.
ARG TESTOPTS="{testopts}"
COPY ./tools/docker/scripts/test_regtest.sh ./
RUN /bin/bash -o pipefail -c " \
    TESTOPTS='${{TESTOPTS}}' \
    ./test_regtest.sh '{arch}' '{version}' |& tee report.log && \
    rm -rf regtesting"
"""
        + print_cached_report()
    )


# ======================================================================================
def regtest_cmake(profile: str, version: str, testopts: str = "") -> str:
    return (
        rf"""
# Install CP2K sources.
WORKDIR /opt/cp2k
COPY ./src ./src
COPY ./data ./data
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .

# Build CP2K with CMake and run regression tests.
ARG TESTOPTS="{testopts}"
COPY ./tools/docker/scripts/build_cp2k_cmake.sh ./tools/docker/scripts/test_regtest_cmake.sh ./
RUN /bin/bash -o pipefail -c " \
    TESTOPTS='${{TESTOPTS}}' \
    ./test_regtest_cmake.sh {profile} {version} |& tee report.log && \
    rm -rf regtesting"
"""
        + print_cached_report()
    )


# ======================================================================================
def regtest_postponed(version: str, arch: str = "local") -> str:
    return (
        install_cp2k(version=version, arch=arch)
        + rf"""
# Postpone running the regression tests until the container is executed.
ARG TESTOPTS
COPY ./tools/docker/scripts/test_regtest.sh ./
ENV TESTOPTS="${{TESTOPTS}}"
CMD ["./test_regtest.sh", "{arch}", "{version}"]

#EOF
"""
    )


# ======================================================================================
def build(version: str, arch: str = "local") -> str:
    return (
        install_cp2k(version=version, arch=arch)
        + rf"""
# Run build test.
COPY ./tools/docker/scripts/test_build.sh .
RUN ./test_build.sh "{arch}" "{version}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def performance(arch: str = "local") -> str:
    return (
        install_cp2k(version="psmp", arch=arch)
        + rf"""
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
        install_cp2k(version=version, arch="local_coverage", revision=True)
        + rf"""
# Run coverage test for {version}.
COPY ./tools/docker/scripts/test_coverage.sh .
RUN ./test_coverage.sh "{version}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def conventions() -> str:
    return (
        rf"""
# Run test for conventions.
WORKDIR /opt/cp2k
COPY ./Makefile .
COPY ./src ./src
COPY ./exts ./exts
COPY ./tools/build_utils ./tools/build_utils
COPY ./tools/conventions ./tools/conventions
COPY ./arch/Linux-x86-64-gfortran.dumpast ./arch/
RUN /bin/bash -ec " \
    ln -vs /opt/cp2k-toolchain/install/arch/local_warn.psmp ./arch/ && \
    source /opt/cp2k-toolchain/install/setup && \
    ./tools/conventions/test_conventions.sh |& tee report.log && \
    rm -rf lib obj exe"
"""
        + print_cached_report()
    )


# ======================================================================================
def manual() -> str:
    return (
        install_cp2k(version="psmp", arch="local", revision=True)
        + rf"""
# Generate manual.
COPY ./docs ./docs
COPY ./tools/input_editing ./tools/input_editing
COPY ./tools/docker/scripts/test_manual.sh .
ARG ADD_EDIT_LINKS=yes
RUN ./test_manual.sh "${{ADD_EDIT_LINKS}}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def precommit() -> str:
    return (
        rf"""
FROM ubuntu:24.04

# Install dependencies.
WORKDIR /opt/cp2k-precommit
COPY ./tools/precommit/ /opt/cp2k-precommit/
RUN ./install_requirements.sh
ENV PATH="/opt/venv/bin:/opt/cp2k-precommit:$PATH"

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
        rf"""
# Install CP2K sources.
WORKDIR /opt/cp2k
COPY ./src ./src
COPY ./data ./data
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .

# Run test for {name}.
COPY ./tools/docker/scripts/build_cp2k_cmake.sh ./tools/docker/scripts/test_{name}.sh ./
RUN ./test_{name}.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def test_without_build(name: str) -> str:
    return (
        rf"""
FROM ubuntu:24.04

# Install dependencies.
WORKDIR /opt/cp2k
COPY ./tools/docker/scripts/install_{name}.sh .
RUN ./install_{name}.sh
ENV PATH="/opt/venv/bin:$PATH"

# Install sources.
ARG GIT_COMMIT_SHA
COPY ./src ./src
COPY ./exts ./exts
COPY ./data ./data
COPY ./docs ./docs
COPY ./tools ./tools
COPY ./tests ./tests
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .
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
ENTRYPOINT []

#EOF
"""


# ======================================================================================
def install_cp2k(
    version: str,
    arch: str,
    revision: bool = False,
    intel: bool = False,
) -> str:
    input_lines = []
    run_lines = []

    if revision:
        input_lines.append("ARG GIT_COMMIT_SHA")
        run_lines.append(
            r'if [ -n "${GIT_COMMIT_SHA}" ] ; then'
            r' echo "git:\${GIT_COMMIT_SHA::7}" > REVISION; fi'
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
        run_lines.append(f"ln -s /opt/cp2k-toolchain /opt/cp2k/tools/toolchain")

    input_block = "\n".join(input_lines)
    run_block = " && \\\n    ".join(run_lines)

    return rf"""
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
def install_dbcsr(profile: str, version: str) -> str:
    return rf"""
# Install DBCSR
COPY ./tools/docker/scripts/install_dbcsr.sh ./
RUN ./install_dbcsr.sh {profile} {version}
"""


# ======================================================================================
def install_deps_toolchain(
    base_image: str = "ubuntu:24.04",
    mpi_mode: str = "mpich",
    with_dbcsr: str = "no",
    with_gcc: str = "system",
    **kwargs: str,
) -> str:
    output = f"\nFROM {base_image}\n\n"
    output += install_toolchain(
        base_image=base_image,
        install_all="",
        mpi_mode=mpi_mode,
        with_dbcsr=with_dbcsr,
        with_gcc=with_gcc,
        **kwargs,
    )
    return output


# ======================================================================================
def install_deps_ubuntu(
    base_image: str = "ubuntu:24.04", gcc_version: int = 13, with_libxsmm: bool = True
) -> str:
    assert gcc_version > 8
    output = rf"""
FROM {base_image}
"""

    if gcc_version > 13:
        output += rf"""
# Add Ubuntu universe repository.
RUN apt-get update -qq && apt-get install -qq --no-install-recommends software-properties-common
RUN add-apt-repository universe
"""

    output += rf"""
# Install Ubuntu packages.
RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
    apt-get update -qq && apt-get install -qq --no-install-recommends \
    cmake \
    less \
    nano \
    make \
    ninja-build \
    wget \
    python3 \
    ca-certificates \
    gcc-{gcc_version} \
    g++-{gcc_version} \
    gfortran-{gcc_version} \
    libfftw3-dev \
    libopenblas-dev \
    libint2-dev \
    libxc-dev \
    libhdf5-dev \
    {"libxsmm-dev" if with_libxsmm else ""} \
    libspglib-f08-dev \
   && rm -rf /var/lib/apt/lists/*

# Create links in /usr/local/bin to overrule links in /usr/bin.
RUN ln -sf /usr/bin/gcc-{gcc_version}      /usr/local/bin/gcc  && \
    ln -sf /usr/bin/g++-{gcc_version}      /usr/local/bin/g++  && \
    ln -sf /usr/bin/gfortran-{gcc_version} /usr/local/bin/gfortran
"""
    return output


# ======================================================================================
def install_deps_ubuntu2004(gcc_version: int = 8) -> str:
    output = rf"""
FROM ubuntu:20.04

# Install Ubuntu packages.
RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
    apt-get update -qq && apt-get install -qq --no-install-recommends \
    cmake \
    gcc-{gcc_version} \
    g++-{gcc_version} \
    gfortran-{gcc_version} \
    libfftw3-dev \
    libopenblas-dev \
    libgsl-dev \
    libhdf5-dev \
   && rm -rf /var/lib/apt/lists/*

# Create links in /usr/local/bin to overrule links in /usr/bin.
RUN ln -sf /usr/bin/gcc-{gcc_version}      /usr/local/bin/gcc  && \
    ln -sf /usr/bin/g++-{gcc_version}      /usr/local/bin/g++  && \
    ln -sf /usr/bin/gfortran-{gcc_version} /usr/local/bin/gfortran

"""
    output += install_toolchain(
        base_image="ubuntu",
        mpi_mode="no",
        with_gcc="system",
        with_cmake="system",
        with_dbcsr="no",
        with_fftw="system",
        with_openblas="system",
        with_gsl="system",
        with_hdf5="system",
        with_libgrpp="no",
        with_libint="install",
        with_libxc="install",
        with_libxsmm="install",
        with_libvori="install",
        with_spglib="no",
    )
    return output


# ======================================================================================
def install_deps_toolchain_intel(
    base_image: str = "intel/hpckit:2024.2.1-0-devel-ubuntu22.04",
    with_ifx: str = "no",
) -> str:
    return rf"""
FROM {base_image}

""" + install_toolchain(
        base_image="ubuntu",
        install_all="",
        with_dbcsr="no",
        with_ifx=with_ifx,
        with_intelmpi="",
        with_mkl="",
        with_libsmeagol="",
        with_libtorch="no",
        with_deepmd="no",
    )


# ======================================================================================
def install_deps_toolchain_nvhpc() -> str:
    return rf"""
FROM ubuntu:22.04

# Install Ubuntu packages.
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    apt-transport-https \
    ca-certificates \
    dirmngr \
    gnupg2 \
    libopenblas-dev \
    make \
    nano \
    python3 \
    wget \
   && rm -rf /var/lib/apt/lists/*

RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK
RUN echo 'deb https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' > /etc/apt/sources.list.d/nvhpc.list

# Install NVIDIA's HPC SDK but only keep the compilers to reduce Docker image size.
RUN apt-get update -qq && \
    apt-get install -qq --no-install-recommends nvhpc-22-11 && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/math_libs && \
    rm -rf /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/comm_libs && \
    rm -rf /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/profilers && \
    rm -rf /opt/nvidia/hpc_sdk/Linux_x86_64/22.11/cuda

ENV PATH ${{PATH}}:/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/bin

# Install CP2K using Linux-x86-64-nvhpc.ssmp.
WORKDIR /opt/cp2k
COPY ./Makefile .
COPY ./src ./src
COPY ./exts ./exts
COPY ./data ./data
COPY ./tests ./tests
COPY ./tools/build_utils ./tools/build_utils
COPY ./tools/regtesting ./tools/regtesting
COPY ./arch/Linux-x86-64-nvhpc.ssmp /opt/cp2k/arch/

# This takes over an hour!
RUN make -j ARCH=Linux-x86-64-nvhpc VERSION=ssmp cp2k
"""


# ======================================================================================
def install_deps_toolchain_cuda(gpu_ver: str) -> str:
    return rf"""
FROM nvidia/cuda:11.8.0-devel-ubuntu22.04

# Setup CUDA environment.
ENV CUDA_PATH /usr/local/cuda
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64

# Disable JIT cache as there seems to be an issue with file locking on overlayfs.
# See also https://github.com/cp2k/cp2k/pull/2337
ENV CUDA_CACHE_DISABLE 1

# Install Ubuntu packages.
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    gfortran                                                          \
    mpich                                                             \
    libmpich-dev                                                      \
   && rm -rf /var/lib/apt/lists/*

""" + install_toolchain(
        base_image="ubuntu",
        mpi_mode="mpich",
        enable_cuda="yes",
        gpu_ver=gpu_ver,
        with_dbcsr="no",
    )


# ======================================================================================
def install_deps_toolchain_hip_cuda(gpu_ver: str) -> str:
    return rf"""
FROM nvidia/cuda:11.8.0-devel-ubuntu22.04

# Setup CUDA environment.
ENV CUDA_PATH /usr/local/cuda
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64
ENV HIP_PLATFORM nvidia
ENV ROCM_VER 4.5.2
ENV HIP_DIR /opt/HIP-rocm-4.5.2
ENV HIPAMD_DIR /opt/hipamd-rocm-4.5.2

# Disable JIT cache as there seems to be an issue with file locking on overlayfs.
# See also https://github.com/cp2k/cp2k/pull/2337
ENV CUDA_CACHE_DISABLE 1

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
    && cd build \
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

# Workaround for HIP installer.
RUN cp -f /opt/hipBLAS-rocm-${{ROCM_VER}}/build/library/src/libhipblas.so /opt/rocm-${{ROCM_VER}}/hipblas/lib/ && \
    cp -f /opt/hipFFT-rocm-${{ROCM_VER}}/build/library/libhipfft.so /opt/rocm-${{ROCM_VER}}/hipfft/lib/

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
ENV LD_LIBRARY_PATH ${{LD_LIBRARY_PATH}}:${{ROCM_PATH}}/lib
ENV HIP_PLATFORM nvidia
RUN hipconfig

""" + install_toolchain(
        base_image="ubuntu",
        mpi_mode="mpich",
        enable_hip="yes",
        gpu_ver=gpu_ver,
        with_dbcsr="no",
    )


# ======================================================================================
def install_deps_toolchain_hip_rocm(gpu_ver: str) -> str:
    return rf"""
FROM rocm/dev-ubuntu-22.04:5.3.2-complete

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
        base_image="ubuntu",
        mpi_mode="mpich",
        enable_hip="yes",
        gpu_ver=gpu_ver,
        with_dbcsr="no",
    )


# ======================================================================================
def install_toolchain(base_image: str, **kwargs: str) -> str:
    install_args = []
    for k, v in kwargs.items():
        k = k.replace("_", "-")
        if v != "":
            install_args.append(f"    --{k}={v} \\")
        else:
            install_args.append(f"    --{k} \\")
    install_args_str = "\n".join(install_args)

    return rf"""
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
# This breaks up the lengthy installation into smaller build steps.
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

COPY ./tools/toolchain/scripts/stage9/ ./scripts/stage9/
RUN  ./scripts/stage9/install_stage9.sh && rm -rf ./build

COPY ./tools/toolchain/scripts/arch_base.tmpl \
     ./tools/toolchain/scripts/generate_arch_files.sh \
     ./scripts/
RUN ./scripts/generate_arch_files.sh && rm -rf ./build
""".lstrip()


# ======================================================================================
def install_deps_spack(version: str) -> str:
    return rf"""
FROM ubuntu:24.04

# Install packages required to build the CP2K dependencies with Spack
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    bzip2 \
    ca-certificates \
    cmake \
    g++ \
    gcc \
    gfortran \
    git \
    gnupg \
    hwloc \
    libhwloc-dev \
    libssh-dev \
    libssl-dev \
    libtool \
    libtool-bin \
    lsb-release \
    make \
    ninja-build \
    patch \
    pkgconf \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    unzip \
    wget \
    xxd \
    xz-utils \
    zstd && rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${{PATH}}"
RUN pip3 install --quiet boto3==1.38.11 google-cloud-storage==3.1.0

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${{NUM_PROCS:-32}}

# Install Spack and Spack packages
WORKDIR /root/spack
ARG SPACK_VERSION
ENV SPACK_VERSION=${{SPACK_VERSION:-1.0.0}}
ARG SPACK_PACKAGES_VERSION
ENV SPACK_PACKAGES_VERSION=${{SPACK_PACKAGES_VERSION:-2025.07.0}}
ARG SPACK_REPO=https://github.com/spack/spack
ENV SPACK_ROOT=/opt/spack-${{SPACK_VERSION}}
ARG SPACK_PACKAGES_REPO=https://github.com/spack/spack-packages
ENV SPACK_PACKAGES_ROOT=/opt/spack-packages-${{SPACK_PACKAGES_VERSION}}
RUN mkdir -p ${{SPACK_ROOT}} \
    && wget -q ${{SPACK_REPO}}/archive/v${{SPACK_VERSION}}.tar.gz \
    && tar -xzf v${{SPACK_VERSION}}.tar.gz -C /opt && rm -f v${{SPACK_VERSION}}.tar.gz \
    && mkdir -p ${{SPACK_PACKAGES_ROOT}} \
    && wget -q ${{SPACK_PACKAGES_REPO}}/archive/v${{SPACK_PACKAGES_VERSION}}.tar.gz \
    && tar -xzf v${{SPACK_PACKAGES_VERSION}}.tar.gz -C /opt && rm -f v${{SPACK_PACKAGES_VERSION}}.tar.gz

ENV PATH="${{SPACK_ROOT}}/bin:${{PATH}}"

# Add Spack packages builtin repository
RUN spack repo add --scope site ${{SPACK_PACKAGES_ROOT}}/repos/spack_repo/builtin

# Find all compilers
RUN spack compiler find

# Find all external packages
RUN spack external find --all --not-buildable

# Add local Spack cache
ARG SPACK_CACHE="s3://spack-cache --s3-endpoint-url=http://localhost:9000"
COPY ./tools/docker/scripts/setup_spack_cache.sh ./
RUN ./setup_spack_cache.sh

# Copy Spack configuration and build recipes
ARG CP2K_VERSION
ENV CP2K_VERSION=${{CP2K_VERSION:-{version}}}
ARG CP2K_BUILD_TYPE
ENV CP2K_BUILD_TYPE=${{CP2K_BUILD_TYPE:-all}}
COPY ./tools/spack/cp2k_deps_${{CP2K_BUILD_TYPE}}_${{CP2K_VERSION}}.yaml ./
COPY ./tools/spack/cp2k_dev_repo ${{SPACK_PACKAGES_ROOT}}/repos/spack_repo/cp2k_dev_repo/
RUN spack repo add --scope site ${{SPACK_PACKAGES_ROOT}}/repos/spack_repo/cp2k_dev_repo/
RUN spack env create myenv cp2k_deps_${{CP2K_BUILD_TYPE}}_${{CP2K_VERSION}}.yaml && \
    spack -e myenv repo list

# Install CP2K dependencies via Spack
RUN spack -e myenv concretize -f
ENV SPACK_ENV_VIEW="${{SPACK_ROOT}}/var/spack/environments/myenv/spack-env/view"
RUN spack -e myenv env depfile -o spack_makefile && \
    make -j${{NUM_PROCS}} --file=spack_makefile SPACK_COLOR=never --output-sync=recurse && \
    cp -ar ${{SPACK_ENV_VIEW}}/bin ${{SPACK_ENV_VIEW}}/include ${{SPACK_ENV_VIEW}}/lib /opt/spack
"""


# ======================================================================================
class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()
        self.content.write(f"#\n")
        self.content.write(f"# This file was created by generate_dockerfiles.py.\n")
        self.content.write(
            f"# Usage: podman build --shm-size=1g -f ./{filename} ../../\n"
        )
        self.content.write(f"#\n")

    def __enter__(self) -> io.StringIO:
        return self.content

    def __exit__(self, exc_type: Any, exc_value: Any, traceback: Any) -> None:
        output_path = Path(__file__).parent / self.filename
        if self.check:
            assert output_path.read_text(encoding="utf8") == self.content.getvalue()
            print(f"File {output_path} is consistent with generator script.")
        else:
            output_path.write_text(self.content.getvalue(), encoding="utf8")
            print(f"Wrote {output_path}")


# ======================================================================================
main()

# EOF
