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
            f.write(install_deps_toolchain(mpi_mode=mpi_mode))
            f.write(regtest("toolchain", version))

    with OutputFile(f"Dockerfile.test_psmp_4ranks", args.check) as f:
        testopts = f"--ompthreads=1 --mpiranks=4"
        f.write(install_deps_toolchain())
        f.write(regtest("toolchain", "psmp", testopts=testopts))

    with OutputFile(f"Dockerfile.test_generic_psmp", args.check) as f:
        f.write(install_deps_toolchain(target_cpu="generic"))
        f.write(regtest("toolchain_generic", "psmp"))

    with OutputFile(f"Dockerfile.test_openmpi-psmp", args.check) as f:
        f.write(install_deps_toolchain(mpi_mode="openmpi"))
        f.write(regtest("toolchain", "psmp"))

    with OutputFile(f"Dockerfile.test_fedora-psmp", args.check) as f:
        f.write(install_deps_toolchain(base_image="fedora:41"))
        f.write(regtest("toolchain", "psmp"))

    for version in "ssmp", "psmp":
        mpi_mode = "intelmpi" if version.startswith("p") else "no"
        f.write(install_deps_toolchain(mpi_mode=mpi_mode))
        with OutputFile(f"Dockerfile.test_intel-ifort-{version}", args.check) as f:
            base_image = "intel/hpckit:2024.2.1-0-devel-ubuntu22.04"
            f.write(install_deps_toolchain_intel(base_image, mpi_mode, with_ifx="no"))
            f.write(regtest("toolchain_intel", version))
        with OutputFile(f"Dockerfile.test_intel-ifx-{version}", args.check) as f:
            base_image = "intel/oneapi-hpckit:2025.2.2-0-devel-ubuntu24.04"
            f.write(install_deps_toolchain_intel(base_image, mpi_mode, with_ifx="yes"))
            f.write(regtest("toolchain_intel", version))

    with OutputFile(f"Dockerfile.test_minimal", args.check) as f:
        f.write(install_deps_ubuntu())
        f.write(regtest("minimal", "ssmp"))

    # Spack/CMake based testers
    with OutputFile(f"Dockerfile.test_spack_psmp", args.check) as f:
        f.write(install_cp2k_spack("psmp", mpi_mode="mpich"))

    for gcc_version in 10, 11, 12, 14:
        with OutputFile(
            f"Dockerfile.test_spack_psmp-gcc{gcc_version}", args.check
        ) as f:
            f.write(
                install_cp2k_spack("psmp", mpi_mode="mpich", gcc_version=gcc_version)
            )

    with OutputFile(f"Dockerfile.test_spack_psmp-4x2", args.check) as f:
        testopts = f"--mpiranks=4 --ompthreads=2"
        f.write(install_cp2k_spack("psmp", mpi_mode="mpich", testopts=testopts))

    with OutputFile(f"Dockerfile.test_spack_openmpi-psmp", args.check) as f:
        f.write(install_cp2k_spack("psmp", mpi_mode="openmpi"))

    with OutputFile(f"Dockerfile.test_spack_ssmp", args.check) as f:
        f.write(install_cp2k_spack("ssmp", mpi_mode="no"))
    # End Spack/CMake based tester

    with OutputFile(f"Dockerfile.test_asan-psmp", args.check) as f:
        f.write(install_deps_toolchain())
        f.write(regtest("toolchain_asan", "psmp"))

    with OutputFile(f"Dockerfile.test_coverage", args.check) as f:
        f.write(install_deps_toolchain())
        f.write(coverage())

    for gcc_version in 8, 9, 10, 11, 12, 13, 14:
        with OutputFile(f"Dockerfile.test_gcc{gcc_version}", args.check) as f:
            # Skip some tests due to bug in LDA_C_PMGB06 functional in libxc <5.2.0.
            testopts = "--skipdir=QS/regtest-rs-dhft" if gcc_version == 8 else ""
            f.write(install_deps_ubuntu(gcc_version=gcc_version))
            f.write(regtest("ubuntu", "ssmp", testopts=testopts))

    with OutputFile("Dockerfile.test_arm64-psmp", args.check) as f:
        base_img = "arm64v8/ubuntu:24.04"
        f.write(install_deps_toolchain(base_img, with_libtorch="no", with_deepmd="no"))
        f.write(regtest("toolchain_arm64", "psmp"))

    with OutputFile(f"Dockerfile.test_performance", args.check) as f:
        f.write(install_deps_toolchain())
        f.write(performance("toolchain"))

    for gpu_ver in "P100", "V100", "A100":
        with OutputFile(f"Dockerfile.test_cuda_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_cuda(gpu_ver=gpu_ver))
            f.write(regtest(f"toolchain_cuda_{gpu_ver}", "psmp"))
        with OutputFile(f"Dockerfile.test_performance_cuda_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_cuda(gpu_ver=gpu_ver))
            f.write(performance(f"toolchain_cuda_{gpu_ver}"))

    for gpu_ver in "Mi50", "Mi100":
        with OutputFile(f"Dockerfile.build_hip_rocm_{gpu_ver}", args.check) as f:
            f.write(install_deps_toolchain_hip_rocm(gpu_ver=gpu_ver))
            f.write(test_build(f"toolchain_hip_{gpu_ver}", "psmp"))

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
            f.write(install_deps_toolchain(mpi_mode="no"))
            f.write(test_3rd_party(name))

    for name in "misc", "doxygen":
        with OutputFile(f"Dockerfile.test_{name}", args.check) as f:
            f.write(test_without_build(name))


# ======================================================================================
def regtest(profile: str, version: str, testopts: str = "") -> str:
    return (
        install_cp2k(profile=profile, version=version)
        + rf"""
# Run regression tests.
ARG TESTOPTS="{testopts}"
COPY ./tests ./tests
COPY ./tools/docker/scripts/test_regtest.sh ./
RUN /bin/bash -o pipefail -c " \
    TESTOPTS='${{TESTOPTS}}' \
    ./test_regtest.sh {profile} {version} |& tee report.log && \
    rm -rf regtesting"
"""
        + print_cached_report()
    )


# ======================================================================================
def test_build(profile: str, version: str) -> str:
    return (
        install_cp2k(profile=profile, version=version)
        + rf"""
# Run build test.
COPY ./tools/docker/scripts/test_build.sh .
RUN ./test_build.sh "{profile}" "{version}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def performance(profile: str) -> str:
    return (
        install_cp2k(profile=profile, version="psmp")
        + rf"""
# Run performance test for {profile}.
COPY ./benchmarks ./benchmarks
COPY ./tools/regtesting ./tools/regtesting
COPY ./tools/docker/scripts/test_performance.sh  \
     ./tools/docker/scripts/plot_performance.py  \
     ./
RUN ./test_performance.sh "{profile}" 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def coverage() -> str:
    return (
        install_cp2k(profile="toolchain_coverage", version="psmp", revision=True)
        + rf"""
# Run coverage test.
COPY ./tests ./tests
COPY ./tools/docker/scripts/test_coverage.sh .
RUN ./test_coverage.sh 2>&1 | tee report.log
"""
        + print_cached_report()
    )


# ======================================================================================
def conventions() -> str:
    return (
        f"""
COPY ./tools/conventions/redirect_gfortran_output.py /usr/bin/
"""
        + install_cp2k(profile="toolchain_conventions", version="psmp")
        + f"""
# Run test for conventions.
COPY ./tools/conventions ./tools/conventions
RUN /bin/bash -ec "./tools/conventions/test_conventions.sh |& tee report.log"
"""
        + print_cached_report()
    )


# ======================================================================================
def manual() -> str:
    return (
        install_cp2k(profile="toolchain", version="psmp", revision=True)
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
        install_cp2k(profile="toolchain", version="ssmp")
        + rf"""
# Run test for {name}.
COPY ./tests ./tests
COPY ./tools/docker/scripts/test_{name}.sh ./
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
COPY ./tools/pao-ml/requirements.txt pao-ml-requirements.txt
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
def install_cp2k(profile: str, version: str, revision: bool = False) -> str:
    output = ""
    if revision:
        output += "\n"
        output += "ARG GIT_COMMIT_SHA\n"
        output += "ENV GIT_COMMIT_SHA=${GIT_COMMIT_SHA}\n"

    output += rf"""
# Install CP2K sources.
WORKDIR /opt/cp2k
COPY ./src ./src
COPY ./data ./data
COPY ./tools/build_utils ./tools/build_utils
COPY ./cmake ./cmake
COPY ./CMakeLists.txt .

# Compile CP2K.
COPY ./tools/docker/scripts/build_cp2k.sh .
RUN ./build_cp2k.sh {profile} {version}
"""
    return output


# ======================================================================================
def install_deps_toolchain(
    base_image: str = "ubuntu:24.04",
    mpi_mode: str = "mpich",
    with_dbcsr: str = "",  # enabled by default
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
def install_deps_ubuntu(gcc_version: int = 13) -> str:
    base_image = "ubuntu:24.04" if gcc_version > 8 else "ubuntu:20.04"
    output = f"\nFROM {base_image}\n"

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
    {"libxsmm-dev" if gcc_version > 8 else ""} \
    {"libspglib-f08-dev" if gcc_version > 8 else ""} \
   && rm -rf /var/lib/apt/lists/*

# Create links in /usr/local/bin to overrule links in /usr/bin.
RUN ln -sf /usr/bin/gcc-{gcc_version}      /usr/local/bin/gcc  && \
    ln -sf /usr/bin/g++-{gcc_version}      /usr/local/bin/g++  && \
    ln -sf /usr/bin/gfortran-{gcc_version} /usr/local/bin/gfortran

# Use toolchain to install DBCSR{"" if gcc_version > 8 else " and CMake"}.
""" + install_toolchain(
        base_image=base_image,
        mpi_mode="no",
        with_dbcsr="",
        with_gcc="system",
        with_cmake="system" if gcc_version > 8 else "",
        with_ninja="system",
        with_openblas="system",
        with_libxc="no",
        with_libint="no",
        with_fftw="no",
        with_libxsmm="no",
        with_spglib="no",
        with_libvori="no",
    )

    return output


# ======================================================================================
def install_deps_toolchain_intel(base_image: str, mpi_mode: str, with_ifx: str) -> str:
    return rf"""
FROM {base_image}

""" + install_toolchain(
        base_image="ubuntu",
        install_all="",
        mpi_mode=mpi_mode,
        with_ifx=with_ifx,
        with_mkl="",
        with_libsmeagol="",
        with_libtorch="no",
        with_deepmd="no",
    )


# ======================================================================================
def install_deps_toolchain_cuda(gpu_ver: str, **kwargs: str) -> str:
    deps = rf"""
FROM nvidia/cuda:12.9.1-devel-ubuntu24.04

# Setup CUDA environment.
ENV CUDA_PATH /usr/local/cuda
ENV LD_LIBRARY_PATH /usr/local/cuda/lib64

# Disable JIT cache as there seems to be an issue with file locking on overlayfs.
# See also https://github.com/cp2k/cp2k/pull/2337
ENV CUDA_CACHE_DISABLE 1

# Install Ubuntu packages.
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    gfortran                                                          \
   && rm -rf /var/lib/apt/lists/*

""" + install_toolchain(
        base_image="ubuntu",
        with_mpich="install",
        mpi_mode="mpich",
        enable_cuda="yes",
        gpu_ver=gpu_ver,
        **kwargs,
    )
    return deps


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

# Remove LTO from Ubuntu's MPICH
RUN sed -i -e 's/-flto=auto//g' -e 's/-ffat-lto-objects//g' \
    /usr/lib/x86_64-linux-gnu/pkgconfig/mpich.pc \
    /usr/bin/*.mpich

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
        with_dbcsr="",
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
     ./tools/toolchain/scripts/generate_cmake_options.sh \
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
""".lstrip()


# ======================================================================================
def install_cp2k_spack(
    version: str,
    mpi_mode: str,
    gcc_version: int = 13,
    testopts: str = "",
) -> str:
    output = rf"""
ARG BASE_IMAGE="ubuntu:24.04"

###### Stage 1: Build CP2K ######

FROM "${{BASE_IMAGE}}" AS build_cp2k

RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    bzip2 \
    ca-certificates \
    cmake \
    g++-{gcc_version} \
    gcc-{gcc_version} \
    gfortran-{gcc_version} \
    git \
    gnupg \
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
    zstd \
    && rm -rf /var/lib/apt/lists/*

# Create links in /usr/local/bin to override tbe links in /usr/bin
RUN ln -sf /usr/bin/gcc-{gcc_version}      /usr/local/bin/gcc && \
    ln -sf /usr/bin/g++-{gcc_version}      /usr/local/bin/g++ && \
    ln -sf /usr/bin/gfortran-{gcc_version} /usr/local/bin/gfortran

ARG SPACK_CACHE="s3://spack-cache --s3-endpoint-url=http://localhost:9000"

# Copy CP2K repository into container
WORKDIR /opt
COPY . cp2k/

# Build CP2K dependencies
WORKDIR /opt/cp2k
RUN /bin/bash -o pipefail -c "source ./make_cp2k.sh -bd_only -cv {version} -mpi {mpi_mode}"

# Build and install CP2K
RUN /bin/bash -o pipefail -c "source ./make_cp2k.sh -cv {version} -mpi {mpi_mode}"

###### Stage 2: Install CP2K ######

FROM "${{BASE_IMAGE}}" AS install_cp2k

# Install required packages
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    g++-{gcc_version} \
    gcc-{gcc_version} \
    gfortran-{gcc_version} \
    python3 \
    && rm -rf /var/lib/apt/lists/*

# Create links in /usr/local/bin to override tbe links in /usr/bin
RUN ln -sf /usr/bin/gcc-{gcc_version}      /usr/local/bin/gcc && \
    ln -sf /usr/bin/g++-{gcc_version}      /usr/local/bin/g++ && \
    ln -sf /usr/bin/gfortran-{gcc_version} /usr/local/bin/gfortran

WORKDIR /opt/cp2k

# Install CP2K dependencies built with spack
COPY --from=build_cp2k /opt/cp2k/spack/spack-1.1.0/opt/spack ./spack/spack-1.1.0/opt/spack

# Install CP2K
COPY --from=build_cp2k /opt/cp2k/install ./install

# Install CP2K regression tests
COPY --from=build_cp2k /opt/cp2k/tests ./tests
COPY --from=build_cp2k /opt/cp2k/src/grid/sample_tasks ./src/grid/sample_tasks

# Install CP2K/Quickstep CI benchmarks
COPY --from=build_cp2k /opt/cp2k/benchmarks/CI ./benchmarks/CI

# Run CP2K regression test
RUN /bin/bash -o pipefail -c "/opt/cp2k/install/bin/entrypoint.sh /opt/cp2k/install/bin/run_tests {testopts}"

# Create entrypoint and finalise container build
WORKDIR /mnt
ENTRYPOINT ["/opt/cp2k/install/bin/entrypoint.sh"]
CMD ["cp2k", "--help"]
"""
    return output


# ======================================================================================
class OutputFile:
    def __init__(self, filename: str, check: bool) -> None:
        self.filename = filename
        self.check = check
        self.content = io.StringIO()
        self.content.write(f"#\n")
        self.content.write(f"# This file was created by generate_dockerfiles.py.\n")
        if "_spack_" in filename or "make_cp2k_" in filename:
            usage = f"./spack_cache_start.sh; podman build --network=host --shm-size=1g -f ./{filename} ../../"
        else:
            usage = f"podman build --shm-size=1g -f ./{filename} ../../"
        self.content.write(f"# Usage: {usage}\n#\n")

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
