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

    for version in "sdbg", "ssmp", "pdbg", "psmp":
        with OutputFile(f"Dockerfile.test_{version}", args.check) as f:
            f.write(toolchain_full() + regtest(version))

        with OutputFile(f"Dockerfile.test_generic_{version}", args.check) as f:
            f.write(toolchain_full(target_cpu="generic") + regtest(version))

    with OutputFile(f"Dockerfile.test_openmpi-psmp", args.check) as f:
        # Also testing --with-gcc=install here, see github.com/cp2k/cp2k/issues/2062 .
        f.write(toolchain_full(mpi_mode="openmpi", with_gcc="install"))
        f.write(regtest("psmp"))

    with OutputFile(f"Dockerfile.test_fedora-psmp", args.check) as f:
        f.write(toolchain_full(base_image="fedora:38") + regtest("psmp"))

    with OutputFile(f"Dockerfile.test_intel-psmp", args.check) as f:
        f.write(toolchain_intel() + regtest("psmp", intel=True))

    with OutputFile(f"Dockerfile.test_nvhpc", args.check) as f:
        f.write(toolchain_nvhpc())

    with OutputFile(f"Dockerfile.test_minimal", args.check) as f:
        f.write(toolchain_full() + regtest("sdbg", "minimal"))

    with OutputFile(f"Dockerfile.test_cmake", args.check) as f:
        f.write(toolchain_full() + regtest_cmake())

    for version in "ssmp", "psmp":
        with OutputFile(f"Dockerfile.test_asan-{version}", args.check) as f:
            f.write(toolchain_full() + regtest(version, "local_asan"))

    for version in "sdbg", "pdbg":
        with OutputFile(f"Dockerfile.test_coverage-{version}", args.check) as f:
            f.write(toolchain_full() + coverage(version))

    for gcc_version in 8, 9, 10, 11, 12:
        with OutputFile(f"Dockerfile.test_gcc{gcc_version}", args.check) as f:
            img = "ubuntu:22.04" if gcc_version > 8 else "ubuntu:20.04"
            f.write(toolchain_ubuntu_nompi(base_image=img, gcc_version=gcc_version))
            # Skip some tests because of bug in LDA_C_PMGB06 functional in libxc <5.2.0.
            f.write(regtest("ssmp", testopts="--skipdir=QS/regtest-rs-dhft"))

    with OutputFile("Dockerfile.test_i386", args.check) as f:
        f.write(toolchain_ubuntu_nompi(base_image="i386/debian:12", libvori=False))
        f.write(regtest("ssmp"))

    with OutputFile("Dockerfile.test_arm64-psmp", args.check) as f:
        f.write(
            toolchain_full(
                base_image="arm64v8/ubuntu:22.04",
                with_libxsmm="no",
                with_libtorch="no",
            )
        )
        f.write(regtest("psmp"))

    with OutputFile(f"Dockerfile.test_performance", args.check) as f:
        f.write(toolchain_full() + performance())

    for gpu_ver in "P100", "V100", "A100":
        with OutputFile(f"Dockerfile.test_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_cuda(gpu_ver=gpu_ver) + regtest("psmp", "local_cuda"))

        with OutputFile(f"Dockerfile.test_hip_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_hip_cuda(gpu_ver=gpu_ver) + regtest("psmp", "local_hip"))

        with OutputFile(f"Dockerfile.test_performance_cuda_{gpu_ver}", args.check) as f:
            f.write(toolchain_cuda(gpu_ver=gpu_ver) + performance("local_cuda"))

    for gpu_ver in "Mi50", "Mi100":
        with OutputFile(f"Dockerfile.test_hip_rocm_{gpu_ver}", args.check) as f:
            # ROCm containers require --device, which is not available for docker build.
            # https://rocmdocs.amd.com/en/latest/ROCm_Virtualization_Containers/ROCm-Virtualization-&-Containers.html#docker-hub
            f.write(toolchain_hip_rocm(gpu_ver=gpu_ver))
            f.write(regtest_postponed("psmp", "local_hip"))

        with OutputFile(f"Dockerfile.build_hip_rocm_{gpu_ver}", args.check) as f:
            f.write(toolchain_hip_rocm(gpu_ver=gpu_ver) + build("psmp", "local_hip"))

    with OutputFile(f"Dockerfile.test_conventions", args.check) as f:
        f.write(toolchain_full() + conventions())

    with OutputFile(f"Dockerfile.test_manual", args.check) as f:
        f.write(toolchain_full() + manual())

    with OutputFile(f"Dockerfile.test_precommit", args.check) as f:
        f.write(precommit())

    for name in "aiida", "ase", "gromacs", "i-pi":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(toolchain_full() + test_3rd_party(name))

    for name in "misc", "doxygen":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(test_without_build(name))

    with OutputFile(f"Dockerfile.gcc_spack", args.check) as f:
        f.write(
            spack_toolchain(
                distro="ubuntu:22.04",
                spec="cp2k@master%gcc build_system=cmake +enable_regtests+cosma+mpi+openmp+sirius+elpa+libxc+libint+plumed+pexsi+spglib ^openblas+fortran ^cosma+scalapack+shared ^dbcsr+mpi+shared",
            )
        )

    with OutputFile(f"Dockerfile.gcc_spack_cuda", args.check) as f:
        f.write(
            spack_toolchain(
                distro="nvidia/cuda:12.2.0-devel-ubuntu22.04",
                spec="cp2k@master%gcc build_system=cmake +libxc+libint+sirius+elpa+plumed+pexsi+spglib+cosma+mpi+openmp+cuda cuda_arch=60  smm=libxsmm ^openblas+fortran ^cosma+scalapack+shared+cuda ^dbcsr+cuda~shared cuda_arch=60 ^sirius+cuda",
            )
        )

    with OutputFile(f"Dockerfile.gcc_spack_rocm", args.check) as f:
        f.write(
            spack_toolchain(
                distro="rocm/dev-ubuntu-22.04:5.6.1-complete",
                spec="cp2k@master%gcc build_system=cmake +sirius +elpa +libxc +libint smm=libxsmm +spglib +cosma +rocm amdgpu_target=gfx90a +pexsi +plumed +openmp ^openblas+fortran ^dbcsr+mpi+rocm~shared+openmp amdgpu_target=gfx906 ^cosma+shared~tests~apps+rocm",
            )
        )


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
def regtest_cmake(testopts: str = "") -> str:
    return (
        rf"""
# Install DBCSR.
WORKDIR /opt/dbcsr
COPY ./tools/docker/scripts/install_dbcsr.sh .
RUN ./install_dbcsr.sh

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
COPY ./tools/docker/scripts/test_regtest_cmake.sh ./
RUN /bin/bash -o pipefail -c " \
    TESTOPTS='${{TESTOPTS}}' \
    ./test_regtest_cmake.sh |& tee report.log && \
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
FROM ubuntu:22.04

# Install dependencies.
WORKDIR /opt/cp2k-precommit
COPY ./tools/precommit/ /opt/cp2k-precommit/
RUN ./install_requirements.sh
ENV PATH="/opt/venv/bin:$PATH"

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
        + rf"""
# Run test for {name}.
COPY ./tools/docker/scripts/test_{name}.sh .
RUN ./test_{name}.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def test_without_build(name: str) -> str:
    return (
        rf"""
FROM ubuntu:22.04

# Install dependencies.
WORKDIR /opt/cp2k
COPY ./tools/docker/scripts/install_{name}.sh .
RUN ./install_{name}.sh

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
def toolchain_full(
    base_image: str = "ubuntu:22.04", with_gcc: str = "system", **kwargs: str
) -> str:
    return f"\nFROM {base_image}\n\n" + install_toolchain(
        base_image=base_image, install_all="", with_gcc=with_gcc, **kwargs
    )


# ======================================================================================
def toolchain_ubuntu_nompi(
    base_image: str = "ubuntu:22.04",
    gcc_version: int = 12,
    libvori: bool = True,
) -> str:
    output = rf"""
FROM {base_image}

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
"""
    if gcc_version > 8:
        output += "    libint2-dev \\\n"
        output += "    libxc-dev \\\n"

    output += rf"""   && rm -rf /var/lib/apt/lists/*

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
        with_fftw="system",
        with_openblas="system",
        with_gsl="system",
        with_hdf5="system",
        with_libint=("system" if gcc_version > 8 else "install"),
        with_libxc=("system" if gcc_version > 8 else "install"),
        with_libxsmm="install",
        with_libvori=("install" if libvori else "no"),
    )
    return output


# ======================================================================================
def toolchain_intel() -> str:
    return rf"""
FROM intel/oneapi-hpckit:2023.2.1-devel-ubuntu22.04

""" + install_toolchain(
        base_image="ubuntu",
        install_all="",
        with_intelmpi="",
        with_mkl="",
        with_libtorch="no",
        with_sirius="no",
    )


# ======================================================================================
def toolchain_nvhpc() -> str:
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
def toolchain_cuda(gpu_ver: str) -> str:
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
        base_image="ubuntu", mpi_mode="mpich", enable_cuda="yes", gpu_ver=gpu_ver
    )


# ======================================================================================
def toolchain_hip_cuda(gpu_ver: str) -> str:
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
        base_image="ubuntu", mpi_mode="mpich", enable_hip="yes", gpu_ver=gpu_ver
    )


# ======================================================================================
def toolchain_hip_rocm(gpu_ver: str) -> str:
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
        base_image="ubuntu", mpi_mode="mpich", enable_hip="yes", gpu_ver=gpu_ver
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


def spack_toolchain(distro: str, spec: str) -> str:
    return rf"""
FROM {distro} as builder
ENV DEBIAN_FRONTEND noninteractive

# only the the two next lines are ubuntu specific
RUN apt-get update -qq
RUN apt-get install -qq --no-install-recommends autoconf autogen automake autotools-dev bzip2 ca-certificates g++ gcc gfortran git less libtool libtool-bin make nano patch pkg-config python3 unzip wget xxd zlib1g-dev cmake gnupg m4 xz-utils libssl-dev libssh-dev openmpi-common libopenmpi-dev 


RUN wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_386 && chmod a+x /usr/local/bin/yq
ENV FORCE_UNSAFE_CONFIGURE 1

ENV PATH="/spack/bin:${{PATH}}"
ENV MPICH_VERSION=4.0.3

# get latest version of spack
RUN git clone https://github.com/spack/spack.git

# set the location of packages built by spack
RUN spack config add config:install_tree:root:/opt/spack

# find all external packages but exclude python
RUN spack external find --all --exclude python

# find compilers
RUN spack compiler find

# tweaking the arguments
RUN yq -i '.compilers[0].compiler.flags.fflags = "-fallow-argument-mismatch"' /root/.spack/linux/compilers.yaml
#RUN spack config add packages:all:target:x86_64
# set amdgpu_target for all packages
# RUN spack config add packages:all:variants:amdgpu_target=gfx90a

# install few dependencies to help with docker caching
RUN spack install cmake@3.27.3
RUN spack install openblas+fortran
RUN spack install libxsmm
RUN spack install libxc
RUN spack install gsl
RUN spack install py-fypp
RUN spack install spglib
RUN spack install libvori

ENV SPEC_OPENBLAS="{spec}"

# instal all dependencies
RUN mkdir /cp2k-src
COPY . /cp2k-src

RUN spack repo add /cp2k-src/tools/spack 
RUN spack env create --with-view /opt/cp2k cp2k-env
RUN spack install --only=dependencies --fail-fast $SPEC_OPENBLAS ^openmpi

# copy source files of the pull request into container
RUN spack --color always -e cp2k-env dev-build -q --source-path /cp2k-src $SPEC_OPENBLAS ^openmpi

# Build CP2K with CMake and run regression tests.
ARG TESTOPTS=""
RUN /bin/bash -o pipefail -c " \
    TESTOPTS='${{TESTOPTS}}' \
    /cp2k-src/tools/docker/scripts/test_regtest_spack.sh |& tee report.log && \
    rm -rf regtesting"

# Output the report if the image is old and was therefore pulled from the build cache.
CMD cat $(find ./report.log -mmin +10) | sed '/^Summary:/ s/$/ (cached)/'
ENTRYPOINT []
"""


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
            print(f"File {output_path} is consistent with generator script.")
        else:
            output_path.write_text(self.content.getvalue(), encoding="utf8")
            print(f"Wrote {output_path}")


# ======================================================================================
main()

# EOF
