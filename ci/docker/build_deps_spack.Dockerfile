# Dockerfile for CP2K continuous integration (CI) runs
#
# A stand-alone build in this folder can be performed with:
# podman build --shm-size=1g -f build_deps_spack.Dockerfile ../../
#
# Author: Matthias Krack (MK)
#
# Stage 1: Create a base image providing the dependencies for building a CP2K binary

ARG BASE_IMAGE="ubuntu:24.04"

FROM "${BASE_IMAGE}" AS build_deps

# Install packages required to build the CP2K dependencies with Spack
RUN apt-get update -qq && apt-get install -qq --no-install-recommends \
    bzip2 \
    ca-certificates \
    cmake \
    g++ gcc gfortran \
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
    zstd \
    && rm -rf /var/lib/apt/lists/*

# Retrieve the number of available CPU cores
ARG NUM_PROCS
ENV NUM_PROCS=${NUM_PROCS:-32}

ARG CP2K_VERSION
ENV CP2K_VERSION=${CP2K_VERSION:-psmp}

# Copy CP2K repository into container
WORKDIR /opt
COPY . cp2k/

# Build CP2K dependencies
WORKDIR /opt/cp2k
RUN /bin/bash -o pipefail -c "source ./make_cp2k.sh -bd_only -cray -cv ${CP2K_VERSION} -dlc -j${NUM_PROCS}"
